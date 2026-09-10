#!/usr/bin/env python3

import sys
from typing import Final

import polars as pl
from lib.common_const import (
    BARCODE,
    HAS_BARCODE,
    HIGH_QUALITY_MAPPED,
    MAPPED,
    MAPPING_QUALITY,
    READ_NAME,
    READS,
    LibraryType,
    Modality,
)
from lib.metrics.metrics import MetricsRecord

# Constants
CHROM: Final[str] = "chrom"
START: Final[str] = "start"
END: Final[str] = "end"
MODAL_START: Final[str] = "modal_start"
MODAL_END: Final[str] = "modal_end"
PLB_READ_NAME: Final[str] = "name"
IS_ANCHORED: Final[str] = "is_anchored"
FLAGS: Final[str] = "flags"

NON_PRIMARY: Final[int] = 0x900
# Thresholds
MAPPING_QUALITY_THRESHOLD: Final[int] = 20

bam_parquet_path = sys.argv[1]  # Bam input path
index_fasta_path = sys.argv[2]  # Genom Fasta path
sample_name = sys.argv[3]
threads = int(sys.argv[4])

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.AMPLICON,
)


full_read_lf = pl.scan_parquet(bam_parquet_path)


positions = (
    full_read_lf.filter(pl.col(CHROM).is_not_null())
    .group_by(CHROM, START)
    .agg(pl.len().alias(READS), pl.mean(MAPPING_QUALITY).alias("mean_mapq"))
)
positions.sink_parquet("positions.parquet")

positions = pl.scan_parquet("positions.parquet")
anchors = (
    positions.sort([READS, START], descending=[True, False])
    .group_by(CHROM, maintain_order=True)
    .agg(pl.col(READS).max(), pl.col(START).first())
    .drop(READS)
)

anchors.sink_parquet("anchors.parquet")
anchors.sink_csv(f"{sample_name}_anchors.csv")

# Modal END anchor per amplicon, the 3'-end analogue of `anchors`. On amplicons
# whose reads sit against a fixed reference END, an indel shifts a read's START by
# the indel length while its END is unchanged, so a START-only filter silently
# drops every indel-bearing read. The END anchor lets those reads back in for
# variant attribution!
end_positions = (
    full_read_lf.filter(pl.col(CHROM).is_not_null())
    .group_by(CHROM, END)
    .agg(pl.len().alias(READS))
)
anchors_end = (
    end_positions.sort([READS, END], descending=[True, False])
    .group_by(CHROM, maintain_order=True)
    .agg(pl.col(END).first())
)

metrics.add("total_reads", full_read_lf.select(pl.len()).collect().item())

# On-target read set, harmonized across ALL amplicon consumers: keep a read if it
# is anchored at the modal START *or* the modal END and clears the mapq threshold.
# An indel keeps one end fixed and shifts the other by its length, so start-OR-end
# accepts indel-bearing reads that a start-only rule silently drops. This single
# definition feeds BOTH cell calling (amplicon_counts) and per-cell variant
# attribution (filtered_amplicon_reads), so every downstream analysis treats
# amplicon reads identically. A read must still anchor at one primer, so genuinely
# off-target reads stay out.
anchored = (
    full_read_lf.filter(pl.col(CHROM).is_not_null())
    .join(anchors.rename({START: MODAL_START}), on=CHROM)
    .join(anchors_end.rename({END: MODAL_END}), on=CHROM)
    .filter((pl.col(FLAGS) & NON_PRIMARY) == 0)
    .filter(pl.col(MAPPING_QUALITY) >= MAPPING_QUALITY_THRESHOLD)
    .filter((pl.col(START) == pl.col(MODAL_START)) | (pl.col(END) == pl.col(MODAL_END)))
)

# Amplicon counts -> cell calling.
filtered = anchored.select(READ_NAME, CHROM, BARCODE)
metrics.add("filtered_reads", filtered.select(pl.len()).collect().item())

counts_lf = (
    filtered.filter(pl.col(BARCODE).is_not_null())
    .group_by(BARCODE, CHROM)
    .agg(pl.len().alias(READS))
)
counts_lf.sink_parquet(f"{Modality.AMPLICON.value}_counts.parquet")

# Read-name allowlist -> per-cell variant attribution (count_variant*.py). Same
# anchored set as the counts above; the two differ only in projected columns.
variant_reads = anchored.select(READ_NAME)
metrics.add("variant_reads", variant_reads.select(pl.len()).collect().item())
variant_reads.sink_parquet("filtered_amplicon_reads.parquet")

metrics.write_records("metrics.parquet")


# Generate upset plot data. IS_ANCHORED uses the same start-OR-end definition as
# the read set above, so the QC classification agrees with what actually feeds
# counts/variants. Left joins keep unmapped / off-target reads (their modal
# start/end come out null -> not anchored, resolved by fill_null(False) below).
raw_joined = (
    full_read_lf.join(anchors.rename({START: MODAL_START}), on=CHROM, how="left")
    .join(anchors_end.rename({END: MODAL_END}), on=CHROM, how="left")
    .with_columns(
        (
            (pl.col(START) == pl.col(MODAL_START)) | (pl.col(END) == pl.col(MODAL_END))
        ).alias(IS_ANCHORED)
    )
)

bool_lf = (
    raw_joined.with_columns(
        pl.when(pl.col(BARCODE).is_not_null())
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(HAS_BARCODE),
        pl.when(pl.col(CHROM).is_not_null())
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(MAPPED),
        pl.when(pl.col(MAPPING_QUALITY) >= MAPPING_QUALITY_THRESHOLD)
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(HIGH_QUALITY_MAPPED),
    )
    .select(IS_ANCHORED, HAS_BARCODE, MAPPED, HIGH_QUALITY_MAPPED)
    .fill_null(False)
    .group_by(pl.all())
    .agg(pl.len().alias(READS))
)

bool_lf.sink_parquet("amplicon_upset_data.parquet")
