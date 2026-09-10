#!/usr/bin/env python3
"""CB-aware per-(cell, amplicon) downsampling keep-list.

Caps each (cell, amplicon) to at most N reads by producing a read-name keep-list.
The caller then subsets the aligned BAM with `samtools view -N`.
Downsampling only ever selects existing reads, so a keep-list is a complete,
lossless description of the subset (without a BAM reconstruction via polarsbio).

Cap strategy: every (cell, amplicon) is capped to min(n, N); nothing is dropped for
being short, so the ADO-prone low-coverage tail is preserved. There is no `>= N` threshold.

The join that defines the eligible read set is the same as per_cell_coverage_qc.py:
read to (amplicon, CB) from the merged reads parquet (CHROM == amplicon contig,
BARCODE == CB), restricted to variant-eligible reads (filtered_amplicon_reads) and to
called cells (top_cells).

Outputs:
  - <sample>_keep_reads.txt: one read_name per line -> `samtools view -N`.
  - <sample>_per_cell_amplicon_subsample_stats.parquet: per (cell, amplicon)
    n_reads_raw / n_reads_used / covered / capped. Once capped, `DP` no longer means
    coverage. The full-depth COUNT_AMPLICON matrix stays the coverage truth.
"""

import sys

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
cap_n = int(sys.argv[5])
seed = int(sys.argv[6])

N_READS_RAW = "n_reads_raw"
N_READS_USED = "n_reads_used"
COVERED = "covered"
CAPPED = "capped"
RANK = "_rank"

keep_list_path = f"{sample_name}_keep_reads.txt"
stats_path = f"{sample_name}_per_cell_amplicon_subsample_stats.parquet"

if cap_n <= 0:
    sys.exit(
        f"cap N (downsample_reads_per_cell_amplicon) must be a positive integer, got {cap_n}"
    )

merged_lf = pl.scan_parquet(merged_reads_path)
filtered_lf = pl.scan_parquet(filtered_amplicon_reads_path)
top_cells_lf = pl.scan_parquet(top_cells_path)

# Distinct (amplicon, CB, read) triples, restricted to variant-eligible reads in
# called cells. read_name is the unit the keep-list caps on, merged_reads holds
# several rows for a read with supplementary alignments, and filtered_amplicon_reads
# is not de-duplicated by the scripts that produce it (count_amplicon.py selects
# read_name off the same multi-row frame).
eligible = (
    merged_lf.filter(pl.col(CHROM).is_not_null())
    .join(filtered_lf.select(READ_NAME), on=READ_NAME, how="inner")
    .join(top_cells_lf.select(BARCODE), on=BARCODE, how="inner")
    .select(CHROM, BARCODE, READ_NAME)
    .unique()
)

# This is for cap-only: Rank the reads within each (amplicon, CB) by a seeded
# per-read hash and keep the cap_n smallest -> a uniform random min(n, cap_n) subset,
# with no threshold (a group of <= cap_n reads keeps all of them). Because the ordering
# key is derived from read_name + seed, the draw is reproducible for a fixed seed and
# independent of the non-deterministic row order polars' scan/join/unique produce,
# so no pre-sort is needed now. read_name is unique within a group, so
# the hashes are distinct there and the rank is unambiguous!
keep = (
    eligible.with_columns(
        pl.col(READ_NAME)
        .hash(seed=seed)
        .rank("ordinal")
        .over(CHROM, BARCODE)
        .alias(RANK)
    )
    .filter(pl.col(RANK) <= cap_n)
    .select(READ_NAME)
)

# read-name list is what `samtools view -N` wants
keep.sink_csv(keep_list_path, include_header=False)

stats = (
    eligible.group_by(CHROM, BARCODE)
    .agg(pl.len().alias(N_READS_RAW))
    .with_columns(
        pl.min_horizontal(pl.col(N_READS_RAW), pl.lit(cap_n, dtype=pl.UInt32)).alias(
            N_READS_USED
        )
    )
    .with_columns(
        (pl.col(N_READS_USED) > 0).alias(COVERED),
        (pl.col(N_READS_RAW) > cap_n).alias(CAPPED),
    )
    .sort(CHROM, BARCODE)
)
stats.sink_parquet(stats_path)

# Keep the run summary from the (per-group, small) stats file for the sake of status
summ = (
    pl.scan_parquet(stats_path)
    .select(
        pl.len().alias("n"),
        pl.col(CAPPED).sum().alias("capped"),
        pl.col(N_READS_RAW).sum().alias("raw"),
        pl.col(N_READS_USED).sum().alias("used"),
    )
    .collect()
    .row(0, named=True)
)
n_groups, n_capped, reads_raw, reads_used = (
    summ["n"],
    summ["capped"],
    summ["raw"],
    summ["used"],
)
print(f"sample\t{sample_name}")
print(f"cap_n\t{cap_n}")
print(f"seed\t{seed}")
print(f"n_cell_amplicons\t{n_groups}")
print(f"n_cell_amplicons_capped\t{n_capped}")
print(f"frac_cell_amplicons_capped\t{(n_capped / n_groups) if n_groups else 0.0:.4f}")
print(f"reads_raw\t{reads_raw}")
print(f"reads_kept\t{reads_used}")
print(f"frac_reads_removed\t{(1 - reads_used / reads_raw) if reads_raw else 0.0:.4f}")
