#!/usr/bin/env python3

import sys
from typing import Final

import polars as pl
from lib.bam_tools.utils import (
    CORRECTED_UMI_TAG,
    GENE_SYMBOL_TAG,
    STAR_BAD_UMI,
    STAR_NO_GENE,
)
from lib.common_const import (
    BARCODE,
    FRACTION,
    HAS_BARCODE,
    HAS_UMI,
    IN_CELL,
    MAPPING_QUALITY,
    PRIMER,
    READ_NAME,
    READS,
    UNMAPPED,
    LibraryType,
    Modality,
)
from lib.metrics.metrics import MetricsRecord
from lib.plotting.amplicon import plot_amplicon_anchors
from lib.plotting.common import export_upset_data_csv, plot_generic_upset_plot

# Constants
CHROM: Final[str] = "chrom"
START: Final[str] = "start"
FLAGS: Final[str] = "flags"
PLB_READ_NAME: Final[str] = "name"

NON_PRIMARY: Final[int] = 0x900
# Thresholds
MAPPING_QUALITY_THRESHOLD: Final[int] = 20

full_read_first_parquet_path = sys.argv[1]  # Bam input path
full_read_second_parquet_path = sys.argv[2]  # Bam input path
top_cells_parquet_path = sys.argv[3]  # Bam input path
index_fasta_path = sys.argv[4]  # Genome Fasta path
sample_name = sys.argv[5]

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.RNA,
    modality=Modality.GENE_EXPRESSION,
)

# Current implementation cannot know which full_read file is from RNA or DNA before loading it.
# TODO: Improve this by adding metadata or standardizing file naming conventions.
temp_lf = pl.scan_parquet(full_read_first_parquet_path)
if CHROM in temp_lf.collect_schema().names():
    full_read_dna_lf = temp_lf
    full_read_rna_lf = pl.scan_parquet(full_read_second_parquet_path)
else:
    full_read_dna_lf = pl.scan_parquet(full_read_second_parquet_path)
    full_read_rna_lf = temp_lf

top_cells = pl.read_parquet(top_cells_parquet_path)

joined = (
    full_read_dna_lf.rename({MAPPING_QUALITY: "mapping_quality_amplicon"})
    .filter((pl.col(FLAGS) & NON_PRIMARY) == 0)
    .join(
        full_read_rna_lf.rename({MAPPING_QUALITY: "mapping_quality_gene_expression"}),
        on=READ_NAME,
        how="inner",
    )
)

positions = (
    full_read_dna_lf.filter(pl.col(CHROM).is_not_null())
    .group_by(CHROM, START, FLAGS)
    .agg(pl.len().alias(READS), pl.mean(MAPPING_QUALITY).alias("mean_mapq"))
)
positions.sink_csv("positions.csv")


positive_filtered_position = positions.filter(
    pl.col(FLAGS) == 0, pl.col(READS) > 10
).collect()
if positive_filtered_position.shape[0] > 2:
    plot_amplicon_anchors(
        positions=positive_filtered_position,
        sample_name=sample_name,
        title="Carryover mapping distributions - Positive Strand",
        prefix="pos",
    )

negative_filtered_position = positions.filter(
    pl.col(FLAGS) == 16, pl.col(READS) > 10
).collect()
if negative_filtered_position.shape[0] > 2:
    plot_amplicon_anchors(
        positions=negative_filtered_position,
        sample_name=sample_name,
        title="Carryover mapping distributions - Negative Strand",
        prefix="neg",
    )

# Generate upset plot data
bool_lf = (
    joined.with_columns(
        pl.when(
            (pl.col(GENE_SYMBOL_TAG) != STAR_NO_GENE)
            & (pl.col(GENE_SYMBOL_TAG).is_not_null())
        )
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("mapping_to_transcriptome"),
        pl.when(pl.col(CORRECTED_UMI_TAG) != STAR_BAD_UMI)
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(HAS_UMI),
        pl.when(pl.col(BARCODE).is_not_null())
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(HAS_BARCODE),
        pl.when(pl.col("mapping_quality_gene_expression") == 255)
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("high_quality_mapped_gene_expression"),
        pl.when(pl.col("mapping_quality_amplicon") >= MAPPING_QUALITY_THRESHOLD)
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("high_quality_mapped_amplicon"),
        pl.when(pl.col(PRIMER) == "GACTTGAGTGGCTGTCGG")
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("has_rna_handle"),
        pl.when(
            (pl.col("mapping_quality_gene_expression") == 0)
            & (
                (pl.col(GENE_SYMBOL_TAG) == STAR_NO_GENE)
                | (pl.col(GENE_SYMBOL_TAG).is_null())
            )
        )
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(UNMAPPED),
        pl.when(
            (pl.col("mapping_quality_gene_expression") == 0)
            & (pl.col(GENE_SYMBOL_TAG) != STAR_NO_GENE)
            & (pl.col(GENE_SYMBOL_TAG).is_not_null())
        )
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("multi-mapped"),
        pl.when(pl.col(BARCODE).is_in(top_cells[BARCODE]))
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias(IN_CELL),
    )
    .select(
        "mapping_to_transcriptome",
        HAS_UMI,
        HAS_BARCODE,
        "high_quality_mapped_gene_expression",
        "high_quality_mapped_amplicon",
        "has_rna_handle",
        UNMAPPED,
        "multi-mapped",
        IN_CELL,
    )
    .fill_null(False)
    .group_by(pl.all())
    .agg(pl.len().alias(READS))
)

plot_generic_upset_plot(
    bool_count=bool_lf.collect(),
    sample_name=sample_name,
    title="Carryover read distribution",
    modularity=Modality.CARRYOVER.value,
    discarded_members=["high_quality_mapped_amplicon", UNMAPPED, "multi-mapped"],
)


# TODO: Add carryover as a modality?
upset_metrics_csv = export_upset_data_csv(
    bool_lf=bool_lf, sample_name=sample_name, prefix="carryover"
)

upset_metrics = pl.read_csv(upset_metrics_csv)
usable_reads = (
    upset_metrics.filter(
        pl.col("has_barcode")
        & pl.col("has_umi")
        & pl.col("mapping_to_transcriptome")
        & pl.col("has_rna_handle")
        & pl.col("high_quality_mapped_gene_expression")
        & ~pl.col("high_quality_mapped_amplicon")
        & ~pl.col("unmapped")
        & ~pl.col("multi-mapped")
    )
    .drop(IN_CELL)
    .sum()
)
usable_reads_in_cell = upset_metrics.filter(
    pl.col("has_barcode")
    & pl.col(IN_CELL)
    & pl.col("has_umi")
    & pl.col("mapping_to_transcriptome")
    & pl.col("has_rna_handle")
    & pl.col("high_quality_mapped_gene_expression")
    & ~pl.col("high_quality_mapped_amplicon")
    & ~pl.col("unmapped")
    & ~pl.col("multi-mapped")
)

fraction_usable_reads = usable_reads.select(FRACTION).item()
total_usable_reads = usable_reads.select(READS).item()
total_usable_reads_in_cell = usable_reads_in_cell.select(READS).item()
fraction_usable_reads_in_cell = usable_reads_in_cell.select(FRACTION).item()

metrics.add("fraction_usable_reads", fraction_usable_reads)
metrics.add("total_usable_reads", total_usable_reads)
metrics.add("fraction_usable_reads_in_cell", fraction_usable_reads_in_cell)
metrics.add("total_usable_reads_in_cell", total_usable_reads_in_cell)
metrics.write_records("metrics.parquet")
