#!/usr/bin/env python3

import sys
from typing import Final

import matplotlib.pyplot as plt
import polars as pl
from lib.bam_tools.utils import (
    CORRECTED_UMI_TAG,
    GENE_SYMBOL_TAG,
    STAR_BAD_UMI,
    STAR_MIN_MAPPING_QUALITY,
    STAR_NO_GENE,
)
from lib.barcode.calling import call_cells
from lib.common_const import (
    BARCODE,
    CONFIDENCE_SCORE,
    COUNTS,
    IS_CELL,
    MAPPING_QUALITY,
    MODALITY,
    N_FEATURES,
    READS,
    TOTAL_COUNTS,
    LibraryType,
    Modality,
)
from lib.dataframe_toolkit import convert_polars_to_adata
from lib.metrics.metrics import MetricsRecord, PerBarcodeMetrics
from lib.pipeline.status import ReasonCode, StatusRecord
from lib.plotting.barcode import (
    experimental_knee_plot,
    plot_feature_count_scatter,
    plot_histograms,
    plot_rank_knee,
)

COUNT_THRESHOLD_DEFAULT_GENE: Final[int] = 200
FEATURE_THRESHOLD_DEFAULT: Final[int] = 1
CONFIDENCE_THRESHOLD_DEFAULT: Final[float] = 0.7

full_read_parquet: str = sys.argv[1]  # Bam input path
sample_name: str = sys.argv[2]
force_cells = sys.argv[3]

if force_cells == "null":
    force_cells = None
else:
    force_cells = int(force_cells)

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.RNA,
    modality=Modality.GENE_EXPRESSION,
)
status = StatusRecord(
    source_id=sample_name,
    library_type=LibraryType.RNA,
    modality=Modality.GENE_EXPRESSION,
)
per_barcode_metrics = PerBarcodeMetrics(
    source_id=sample_name,
    library_type=LibraryType.RNA,
    modality=Modality.GENE_EXPRESSION,
)

full_read_lf = pl.scan_parquet(full_read_parquet)
total_reads = full_read_lf.select(pl.len()).collect().item()


valid_reads_filters = (
    pl.col(GENE_SYMBOL_TAG) != STAR_NO_GENE,
    pl.col(CORRECTED_UMI_TAG) != STAR_BAD_UMI,
    pl.col(MAPPING_QUALITY) == STAR_MIN_MAPPING_QUALITY,
    pl.col(BARCODE).is_not_null(),
)

feature_filtered_counts_lf = (
    full_read_lf.filter(*valid_reads_filters)
    .group_by(BARCODE, GENE_SYMBOL_TAG, CORRECTED_UMI_TAG)
    .agg(pl.len().alias(READS))
    .group_by(BARCODE, GENE_SYMBOL_TAG)
    .agg(pl.len().alias(COUNTS))
)

sorted_per_cell_total = (
    feature_filtered_counts_lf.group_by(BARCODE)
    .agg(pl.sum(COUNTS))
    .sort(COUNTS, descending=True)
).collect(engine="streaming")


barcode_metrics = (
    feature_filtered_counts_lf.group_by(BARCODE)
    .agg(
        pl.col(GENE_SYMBOL_TAG).n_unique().alias(N_FEATURES),
        pl.sum(COUNTS).alias(TOTAL_COUNTS),
    )
    .collect(engine="streaming")
)

per_barcode_metrics.join_metrics(
    barcode_metrics.select(BARCODE, N_FEATURES, TOTAL_COUNTS)
)

barcode_scores = call_cells(
    barcode_metrics=barcode_metrics,
    count_floor=COUNT_THRESHOLD_DEFAULT_GENE,
    feature_count_floor=FEATURE_THRESHOLD_DEFAULT,
    min_confidence=CONFIDENCE_THRESHOLD_DEFAULT,
)

per_barcode_metrics.join_metrics(
    barcode_scores.select(BARCODE, IS_CELL, CONFIDENCE_SCORE)
)
confidence_summary = barcode_scores.group_by(IS_CELL).agg(
    pl.mean(CONFIDENCE_SCORE).alias("mean_confidence"),
    pl.median(CONFIDENCE_SCORE).alias("median_confidence"),
)


filtered_barcodes_mean_confidence = (
    confidence_summary.filter(pl.col(IS_CELL))
    .select(
        pl.col("mean_confidence").first()
    )  # Trick to get a none from polars if the filtered barcodes are empty
    .item()
)
filtered_barcodes_median_confidence = (
    confidence_summary.filter(pl.col(IS_CELL))
    .select(pl.col("median_confidence").first())
    .item()
)

metrics.add("filtered_barcodes_mean_confidence", filtered_barcodes_mean_confidence)
metrics.add("filtered_barcodes_median_confidence", filtered_barcodes_median_confidence)

experimental_knee_plot(
    data=barcode_scores,
    sample_name=sample_name,
    ranking_col=TOTAL_COUNTS,
    modality=Modality.GENE_EXPRESSION.value,
)

plot_feature_count_scatter(
    data=barcode_scores,
    sample_name=sample_name,
    confidence_threshold=CONFIDENCE_THRESHOLD_DEFAULT,
)

# Use force cells from input samplesheet
if force_cells is not None:
    top_cells_stats = barcode_scores.sort(TOTAL_COUNTS, descending=True).head(
        force_cells
    )
else:
    top_cells_stats = barcode_scores.filter(pl.col(IS_CELL))

top_cells = top_cells_stats.select(BARCODE)

over_1000_cells = sorted_per_cell_total.filter(pl.col(COUNTS) > 1000).select(BARCODE)
metrics.add("total_cells_over_1000_umis", over_1000_cells.shape[0])


metrics.add("total_raw_reads", total_reads)

total_reads_in_cells = (
    full_read_lf.join(top_cells.lazy(), on=BARCODE).select(pl.len()).collect().item()
)
metrics.add("total_reads_in_cells", total_reads_in_cells)
filtered_cells = top_cells.shape[0]
metrics.add("filtered_cells", filtered_cells)
fraction_reads_in_cells = total_reads_in_cells / total_reads if total_reads > 0 else 0
metrics.add("fraction_reads_in_cells", fraction_reads_in_cells)

feature_filtered_counts_lf.sink_parquet("raw_counts.parquet")

feature_filtered_counts_lf.join(top_cells.lazy(), on=BARCODE).sink_parquet(
    "filtered_reads.parquet"
)
if filtered_cells == 0:
    status.add(
        ReasonCode.NO_CELLS,
        "No cells were found after filtering.",
    )
else:
    top_cells.write_parquet(f"{sample_name}_top_cells.parquet")


adata = convert_polars_to_adata(
    feature_filtered_counts_lf.collect(),
    barcode_index_name=BARCODE,
    feature_index_name=GENE_SYMBOL_TAG,
    values=COUNTS,
    annotations=barcode_scores.select(BARCODE, CONFIDENCE_SCORE, IS_CELL)
    .to_pandas()
    .set_index(BARCODE),
)
adata.uns[MODALITY] = Modality.GENE_EXPRESSION.value


adata.write_h5ad("raw_gene_counts.h5ad")
metrics.add(
    "total_barcodes",
    full_read_lf.select(BARCODE).unique().select(pl.len()).collect().item(),
)


plt.figure(figsize=(6, 4), dpi=120)
umi_threshold = plot_rank_knee(
    sorted_per_cell_total,
    knee_ncells=top_cells.shape[0],
    sample_name=sample_name,
    value_type=COUNTS,
)

plot_histograms(
    sorted_per_cell_total,
    sample_name=sample_name,
    umi_threshold=umi_threshold,
    data_type=COUNTS,
)


metrics.write_records("metrics.parquet")
per_barcode_metrics.write_records("per_barcode_metrics.parquet")
status.write_records("status.parquet")
