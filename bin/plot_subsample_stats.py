#!/usr/bin/env python3
"""Before/after picture of the per-(cell, amplicon) read cap for one sample.

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

stats_path = sys.argv[1]
sample_name = sys.argv[2]
cap_n = int(sys.argv[3])

N_READS_RAW = "n_reads_raw"
N_READS_USED = "n_reads_used"
CAPPED = "capped"

OVERFLOW_PCTL = 99.5
OVERFLOW_MIN_CAP_MULT = 2
TICK_STEPS = (10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000)
MAX_TICKS = 5

N_BINS = 40
BAR, CAPC, INK, MUTED, FAINT = "#5b8db8", "#d63b3b", "#1e293b", "#64748b", "#94a3b8"

out_png = f"{sample_name}_subsample_before_after.png"

if cap_n <= 0:
    sys.exit(f"cap N must be a positive integer, got {cap_n}")

df = pl.read_parquet(stats_path)
raw = df[N_READS_RAW].to_numpy().astype(float)
used = df[N_READS_USED].to_numpy().astype(float)
n_capped = int(df[CAPPED].sum()) if df.height else 0

print(
    f"sample\t{sample_name}\ncap_n\t{cap_n}\npairs\t{raw.size}\n"
    f"reads_raw\t{int(raw.sum())}\nreads_kept\t{int(used.sum())}\n"
    f"capped_pairs\t{n_capped}"
)

if raw.size == 0:
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    ax.text(0.5, 0.5, "no (cell, amplicon) pairs", ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(out_png)
    sys.exit(0)

target = max(OVERFLOW_MIN_CAP_MULT * cap_n, np.percentile(raw, OVERFLOW_PCTL))
step = next(
    (s for s in TICK_STEPS if target / s <= MAX_TICKS),
    math.ceil(target / MAX_TICKS),
)
xmax = int(step * math.ceil(target / step))
ticks = list(range(0, xmax + 1, int(step)))
bins = np.linspace(0, xmax, N_BINS + 1)

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3), sharey=True, dpi=120)

for ax, data, title in [
    (axes[0], raw, "as sequenced"),
    (axes[1], used, "after capping"),
]:
    ax.hist(
        np.clip(data, None, xmax - 1),
        bins=bins,
        color=BAR,
        edgecolor="white",
        linewidth=0.6,
    )
    ax.axvline(cap_n, color=CAPC, lw=1.8, ls="--")
    ax.set_title(title, fontsize=14, fontweight="600", color=INK, pad=12)
    ax.set_xlabel("reads for one amplicon in one cell", fontsize=11.5, color=MUTED)
    ax.set_xlim(0, xmax)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks[:-1]] + [f"{xmax}+"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#cbd5e1")
    ax.tick_params(colors=MUTED, labelsize=10)

top = axes[0].get_ylim()[1]

axes[0].text(
    cap_n + xmax * 0.02,
    top * 0.94,
    f"cap = {cap_n}",
    color=CAPC,
    fontsize=11,
    fontweight="600",
    va="top",
)
axes[0].set_ylabel("number of cell-amplicon pairs", fontsize=11.5, color=MUTED)
axes[0].text(
    xmax * 0.34,
    top * 0.55,
    f"{100 * n_capped / raw.size:.1f}% of pairs sit above the cap\n"
    f"capping drops {100 * (1 - used.sum() / raw.sum()):.1f}% of reads",
    fontsize=11,
    color=MUTED,
    ha="left",
)
axes[0].text(
    xmax * 0.34,
    top * 0.36,
    f"last bar is everything above {xmax}; the busiest is {int(raw.max())}",
    fontsize=9.5,
    color=FAINT,
    ha="left",
)
axes[1].text(
    cap_n + xmax * 0.06,
    top * 0.55,
    "everything above the cap\ncollapses to the same depth",
    fontsize=11,
    color=MUTED,
    ha="left",
)

fig.suptitle(sample_name, fontsize=11, color=MUTED, y=1.02)
fig.tight_layout(pad=1.4)
fig.savefig(out_png, bbox_inches="tight", facecolor="white", pad_inches=0.2)
print("wrote", out_png)
