#!/usr/bin/env python3
"""Per-(cell, amplicon) read depth as sequenced, with the variant-calling cap, for one sample.

In:  <sample>_per_cell_amplicon_subsample_stats.parquet (chrom, barcode,
     n_reads_raw, n_reads_used, covered, capped), sample name, cap N.
Out: <sample>_subsample_before_after.png -> qc

Run: PYTHONPATH=bin python3 bin/plot_subsample_stats.py <stats.parquet> <sample> 100
"""

import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.ticker import FuncFormatter

stats_path = sys.argv[1]
sample_name = sys.argv[2]
cap_n = int(sys.argv[3])

N_READS_RAW = "n_reads_raw"
N_READS_USED = "n_reads_used"

BINS_PER_DECADE = 10
HEADROOM = 1.15
BAR, CAPC, INK, MUTED = "#5b8db8", "#d63b3b", "#1e293b", "#64748b"

out_png = f"{sample_name}_subsample_before_after.png"

if cap_n <= 0:
    sys.exit(f"cap N must be a positive integer, got {cap_n}")

df = pl.read_parquet(stats_path)
raw = df[N_READS_RAW].to_numpy()
used = df[N_READS_USED].to_numpy()
n_at_cap = int((raw >= cap_n).sum())

print(
    f"sample\t{sample_name}\ncap_n\t{cap_n}\npairs\t{raw.size}\n"
    f"reads_raw\t{int(raw.sum())}\nreads_kept\t{int(used.sum())}\n"
    f"pairs_at_cap\t{n_at_cap}"
)

if raw.size == 0:
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    ax.text(0.5, 0.5, "no (cell, amplicon) pairs", ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(out_png)
    sys.exit(0)

# log-spaced
k = np.arange(
    -math.ceil(math.log10(cap_n) * BINS_PER_DECADE),
    math.ceil(math.log10(max(raw.max(), 2 * cap_n) / cap_n) * BINS_PER_DECADE) + 1,
)
edges = np.unique(np.round(cap_n * 10.0 ** (k / BINS_PER_DECADE)).astype(int))
edges = edges[edges >= 1]
counts, _ = np.histogram(raw, bins=edges)

fig, ax = plt.subplots(figsize=(8.4, 4.3), dpi=120)
ax.bar(
    edges[:-1],
    counts,
    width=np.diff(edges),
    align="edge",
    color=BAR,
    edgecolor="white",
    linewidth=0.6,
)
ax.axvline(cap_n, color=CAPC, lw=1.8, ls="--")
ax.set_xscale("log")
ax.set_xlim(1, edges[-1])
ax.set_ylim(0, max(counts.max(), 1) * HEADROOM)
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.set_xlabel(
    "reads for one amplicon in one cell, as sequenced (log scale)",
    fontsize=11.5,
    color=MUTED,
)
ax.set_ylabel("number of cell-amplicon pairs", fontsize=11.5, color=MUTED)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color("#cbd5e1")
ax.tick_params(colors=MUTED, labelsize=10)

ax.text(
    cap_n * 1.08,
    ax.get_ylim()[1] * 0.97,
    f"cap = {cap_n}",
    color=CAPC,
    fontsize=11,
    fontweight="600",
    va="top",
)

# the stats only hold pairs with >= 1 read, so a pair with no reads is in no denominator
ax.set_title(
    f"{100 * n_at_cap / raw.size:.1f}% of cell-amplicon pairs with reads reach the cap "
    f"of {cap_n}.\nCapping keeps {100 * used.sum() / raw.sum():.1f}% of reads for "
    "variant calling.",
    fontsize=11,
    color=MUTED,
    loc="left",
)

fig.suptitle(sample_name, fontsize=12, color=INK, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(out_png, bbox_inches="tight", facecolor="white", pad_inches=0.2)
print("wrote", out_png)
