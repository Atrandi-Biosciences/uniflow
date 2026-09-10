#!/usr/bin/env python3
"""Per-(cell, amplicon) coverage QC

Generates the per-cell coverage picture so a full-depth run can answer whether
to downsample and at what cap N. This runs on every amplicon run regardless of
whether per-cell/amplicon downsampling is enabled.

Inputs mirror the CB-aware keep-list join in subsample_keeplist.py, that is
read->(amplicon, CB) comes from the merged reads parquet (CHROM == amplicon
contig, BARCODE == CB). The read set is restricted to the variant-eligible reads
(filtered_amplicon_reads, read_name only) and to called cells (top_cells).
Nothing is downsampled here; this only measures what is seen and what a cap would
do.

Outputs (all derived from the one per-(cell, amplicon) depth parquet):
  * <sample>_per_cell_amplicon_depth.parquet     per-(cell, amplicon) read counts   -> dev
  * <sample>_per_cell_amplicon_depth_summary.txt  per-amplicon frac_below_cap etc.  -> dev
  * <sample>_cap_impact_sweep.csv                 % capped / % reads removed vs N    -> dev
  * <sample>_depth_histogram.csv                  read-count bin distribution        -> dev
  * <sample>_per_cell_amplicon_rank_knee_pooled.png    rank-knee, all cell-amplicons -> qc
  * <sample>_per_cell_amplicon_rank_knee_faceted.png   rank-knee, per-amplicon facet -> qc
  * <sample>_per_cell_total_rank_knee.png              rank-knee, per-cell total      -> qc
  * <sample>_cap_impact_curve.png                      % reads kept / % uncapped vs N -> qc
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotnine as gg
import polars as pl
from lib.common_const import (
    BARCODE,
    CHROM,
    READ_NAME,
)

# full_read_*.parquet: read_name, chrom, barcode
merged_reads_path = sys.argv[1]
# filtered_amplicon_reads.parquet: read_name only
filtered_amplicon_reads_path = sys.argv[2]
# *_top_cells.parquet: barcode
top_cells_path = sys.argv[3]
sample_name = sys.argv[4]
# currently-configured N (marker on plots/tables)
cap_n = int(sys.argv[5])

N_READS = "n_reads"
RANK = "rank"

# Candidate downsampling caps for the sweep table + cap-impact curve.
CANDIDATE_GRID = (
    [1, 2, 5, 10, 20, 30, 40] + list(range(50, 1001, 50)) + list(range(1100, 3001, 100))
)

# Coarser subset of CANDIDATE_GRID drawn as reference lines on the pooled rank-knee
# (the full grid is for the table/curve, not the plot).
KNEE_REF_LINES = [1, 500, 1000, 1500, 2000, 2500, 3000]

# Read-count bins for the depth histogram
HIST_BINS = [
    ("1-5", 1, 5),
    ("6-10", 6, 10),
    ("11-20", 11, 20),
    ("21-50", 21, 50),
    ("51-100", 51, 100),
    ("101-200", 101, 200),
    ("200+", 201, None),
]


# ---------------------------------------------------------------------------
# per-(cell, amplicon) depth
# ---------------------------------------------------------------------------
merged_lf = pl.scan_parquet(merged_reads_path)
filtered_lf = pl.scan_parquet(filtered_amplicon_reads_path)
top_cells_lf = pl.scan_parquet(top_cells_path)

# TODO: unify this filter with the other scripts use this such as count_amplicon* scripts
# read -> (amplicon, CB), restricted to variant-eligible reads in called cells.

per_cell_amplicon = (
    merged_lf.filter(pl.col(CHROM).is_not_null())
    .join(filtered_lf.select(READ_NAME), on=READ_NAME, how="inner")
    .join(top_cells_lf.select(BARCODE), on=BARCODE, how="inner")
    .group_by(CHROM, BARCODE)
    .agg(pl.col(READ_NAME).n_unique().alias(N_READS))
)

dist = per_cell_amplicon.collect()
dist.write_parquet(f"{sample_name}_per_cell_amplicon_depth.parquet")

n_groups = dist.height
# int64: n_unique() is unsigned, and `counts - n` underflowed...
counts = dist.get_column(N_READS).to_numpy().astype(np.int64)


# ---------------------------------------------------------------------------
# summary (per-amplicon breakdown)
# ---------------------------------------------------------------------------
if n_groups == 0:
    summary_lines = [
        f"sample\t{sample_name}",
        f"cap_n\t{cap_n}",
        "n_cell_amplicons\t0",
        "WARNING: no (cell, amplicon) groups found; check inputs.",
    ]
else:
    n_below = int((counts < cap_n).sum())
    frac_below = n_below / n_groups
    summary_lines = [
        f"sample\t{sample_name}",
        f"cap_n\t{cap_n}",
        f"n_cell_amplicons\t{n_groups}",
        f"n_reads_total\t{int(counts.sum())}",
        f"reads_min\t{int(counts.min())}",
        f"reads_p50\t{float(np.median(counts)):.1f}",
        f"reads_p90\t{float(np.quantile(counts, 0.90)):.1f}",
        f"reads_max\t{int(counts.max())}",
        f"n_below_cap\t{n_below}",
        f"frac_below_cap\t{frac_below:.4f}",
    ]

    # Per-amplicon median so a single flat N can be judged against amplicons that
    # differ by orders of magnitude in efficiency.
    per_amplicon = (
        dist.group_by(CHROM)
        .agg(
            pl.len().alias("n_cells"),
            pl.col(N_READS).median().alias("median_reads"),
            pl.col(N_READS).max().alias("max_reads"),
            (pl.col(N_READS) < cap_n).mean().alias("frac_below_cap"),
        )
        .sort(CHROM)
    )
    summary_lines.append("")
    summary_lines.append(
        "per_amplicon\tn_cells\tmedian_reads\tmax_reads\tfrac_below_cap"
    )
    for row in per_amplicon.iter_rows(named=True):
        summary_lines.append(
            f"{row[CHROM]}\t{row['n_cells']}\t{row['median_reads']:.1f}\t"
            f"{row['max_reads']}\t{row['frac_below_cap']:.4f}"
        )

summary = "\n".join(summary_lines) + "\n"
with open(f"{sample_name}_per_cell_amplicon_depth_summary.txt", "w") as fh:
    fh.write(summary)
print(summary)


# ---------------------------------------------------------------------------
# table 1: cell cap-impact sweep  (respecting observed max set above)
# ---------------------------------------------------------------------------
grid = [] if n_groups == 0 else list(CANDIDATE_GRID)

total_reads = int(counts.sum()) if n_groups else 0
sweep_rows = []
for n in grid:
    n_capped = int((counts > n).sum())
    reads_removed = int(np.clip(counts - n, 0, None).sum())
    sweep_rows.append(
        {
            "N": n,
            "n_cell_amplicons_capped": n_capped,
            "pct_cell_amplicons_capped": round(100.0 * n_capped / n_groups, 4),
            "reads_total": total_reads,
            "reads_removed": reads_removed,
            "pct_reads_removed": (
                round(100.0 * reads_removed / total_reads, 4) if total_reads else 0.0
            ),
        }
    )
sweep_schema = {
    "N": pl.Int64,
    "n_cell_amplicons_capped": pl.Int64,
    "pct_cell_amplicons_capped": pl.Float64,
    "reads_total": pl.Int64,
    "reads_removed": pl.Int64,
    "pct_reads_removed": pl.Float64,
}
sweep_df = pl.DataFrame(sweep_rows, schema=sweep_schema)
sweep_df.write_csv(f"{sample_name}_cap_impact_sweep.csv")


# ---------------------------------------------------------------------------
# table 2: depth histogram  (fixed read-count bins)
# ---------------------------------------------------------------------------
hist_rows = []
cumulative = 0
for label, lo, hi in HIST_BINS:
    if n_groups == 0:
        in_bin = 0
    elif hi is None:
        in_bin = int((counts >= lo).sum())
    else:
        in_bin = int(((counts >= lo) & (counts <= hi)).sum())
    cumulative += in_bin
    pct = round(100.0 * in_bin / n_groups, 4) if n_groups else 0.0
    cum_pct = round(100.0 * cumulative / n_groups, 4) if n_groups else 0.0
    hist_rows.append(
        {
            "bin": label,
            "n_cell_amplicons": in_bin,
            "pct": pct,
            "cumulative_pct": cum_pct,
        }
    )
hist_df = pl.DataFrame(
    hist_rows,
    schema={
        "bin": pl.Utf8,
        "n_cell_amplicons": pl.Int64,
        "pct": pl.Float64,
        "cumulative_pct": pl.Float64,
    },
)
hist_df.write_csv(f"{sample_name}_depth_histogram.csv")


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def _placeholder(out_png: str, msg: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    ax.text(0.5, 0.5, msg, ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(out_png)
    plt.close(fig)


def _rank_knee_mpl(values, title, out_png, ref_lines=None, highlight=None):
    """Descending rank (x, log) vs read count (y, log)"""
    y = np.sort(np.asarray(values))[::-1]
    fig, ax = plt.subplots(figsize=(8, 6), dpi=120)
    ranks = np.arange(1, len(y) + 1)
    ax.plot(ranks, y, lw=1.2, color="#1f77b4")
    for ref in ref_lines or []:
        ax.axhline(ref, color="0.75", ls=":", lw=0.8)
    if highlight is not None and highlight > 0:
        ax.axhline(
            highlight,
            color="crimson",
            ls="--",
            lw=1.2,
            label=f"configured N = {highlight}",
        )
        ax.legend(frameon=False)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Rank (log)")
    ax.set_ylabel("Reads (log)")
    ax.set_title(f"{title}\n{sample_name}")
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


pooled_png = f"{sample_name}_per_cell_amplicon_rank_knee_pooled.png"
faceted_png = f"{sample_name}_per_cell_amplicon_rank_knee_faceted.png"
percell_png = f"{sample_name}_per_cell_total_rank_knee.png"
curve_png = f"{sample_name}_cap_impact_curve.png"

if n_groups == 0:
    for png in (pooled_png, faceted_png, percell_png, curve_png):
        _placeholder(png, "no (cell, amplicon) groups")
else:
    # 1. pooled rank-knee over all (cell, amplicon) groups. Reference lines are the
    # coarse KNEE_REF_LINES subset (clipped to max), not the full sweep grid.
    _rank_knee_mpl(
        counts,
        "Per-(cell, amplicon) coverage rank curve",
        pooled_png,
        ref_lines=[n for n in KNEE_REF_LINES if n <= int(counts.max())],
        highlight=cap_n,
    )

    # 2. faceted rank-knee, one panel per amplicon (assumes <= 60 amplicons for now).
    ranked = dist.with_columns(
        pl.col(N_READS).rank("ordinal", descending=True).over(CHROM).alias(RANK)
    )
    n_amplicons = ranked.get_column(CHROM).n_unique()
    ncol = min(6, max(1, n_amplicons))
    nrow = int(np.ceil(n_amplicons / ncol))
    facet_plot = (
        gg.ggplot(ranked, gg.aes(x=RANK, y=N_READS))
        + gg.geom_point(size=0.5, alpha=0.6, color="#1f77b4")
        + gg.geom_hline(yintercept=cap_n, color="crimson", linetype="dashed")
        + gg.scale_x_log10()
        + gg.scale_y_log10()
        + gg.facet_wrap(CHROM, scales="free_x", ncol=ncol)
        + gg.labs(
            x="Rank within amplicon (log)",
            y="Reads (log)",
            title=f"Per-cell coverage per amplicon (dashed = configured N={cap_n})",
            subtitle=sample_name,
        )
        + gg.theme(figure_size=(3.2 * ncol, 2.6 * nrow))
    )
    facet_plot.save(faceted_png, dpi=120, verbose=False)

    # 3. per-cell total rank-knee (reads summed over amplicons per cell).
    per_cell_total = dist.group_by(BARCODE).agg(pl.col(N_READS).sum().alias(N_READS))
    _rank_knee_mpl(
        per_cell_total.get_column(N_READS).to_numpy(),
        "Per-cell total coverage rank curve",
        percell_png,
    )

    # 4. cell cap-impact curve: % reads kept and % cell-amplicons
    # kept whole (uncapped) vs N
    fig, ax = plt.subplots(figsize=(8, 6), dpi=120)
    ns = sweep_df.get_column("N").to_numpy()
    pct_reads_kept = 100.0 - sweep_df.get_column("pct_reads_removed").to_numpy()
    pct_uncapped = 100.0 - sweep_df.get_column("pct_cell_amplicons_capped").to_numpy()
    ax.plot(ns, pct_reads_kept, "-o", color="#2ca02c", label="% reads kept")
    ax.plot(
        ns,
        pct_uncapped,
        "-s",
        color="#1f77b4",
        label="% cell-amplicons kept whole",
    )
    ax.axvline(cap_n, color="crimson", ls="--", lw=1.0, label=f"configured N = {cap_n}")
    ax.set_xscale("log")
    ax.set_xlabel("Cap N (log)")
    ax.set_ylabel("Percent kept")
    ax.set_ylim(0, 100)
    ax.set_title(f"Cap impact vs N\n{sample_name}")
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(curve_png)
    plt.close(fig)
