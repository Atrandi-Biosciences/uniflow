#!/usr/bin/env python3


import gc
import os
import sys
from typing import Final

import polars as pl
from lib.bam_tools.utils import (
    CORRECTED_UMI_TAG,
    GENE_SYMBOL_TAG,
    STAR_BAD_UMI,
    STAR_NO_GENE,
)
from lib.common_const import BARCODE, MAPPING_QUALITY, LibraryType, Modality
from lib.library_qc import (
    REFERENCE_DATASETS,
    STAR_MIN_MAPPING_QUALITY,
    get_saturation_curve_cells_reverse,
    plot_comparative_saturation_curves,
)
from lib.metrics.metrics import MetricsRecord

jemalloc_conf = "dirty_decay_ms:500,muzzy_decay_ms:-1"
if os.environ.get("POLARS_THP") == "1":
    jemalloc_conf += ",thp:always,metadata_thp:always"
if override := os.environ.get("_RJEM_MALLOC_CONF"):
    jemalloc_conf += "," + override
os.environ["_RJEM_MALLOC_CONF"] = jemalloc_conf

full_read_parquet: str = sys.argv[1]  # Bam input path
top_cells_path: str = sys.argv[2]
sample_name: str = sys.argv[3]

# Constants
FILTERED_INPUT_DATA: Final[str] = "filtered"

# Load data
top_cells_df = pl.read_parquet(top_cells_path)
full_read_df = pl.scan_parquet(full_read_parquet)

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.RNA,
    modality=Modality.GENE_EXPRESSION,
)

filtered = full_read_df.join(top_cells_df.lazy(), on=BARCODE).filter(
    pl.col(MAPPING_QUALITY) == STAR_MIN_MAPPING_QUALITY,
    pl.col(GENE_SYMBOL_TAG) != STAR_NO_GENE,
    pl.col(CORRECTED_UMI_TAG) != STAR_BAD_UMI,
)

del full_read_df
gc.collect()

filtered_df_metrics = get_saturation_curve_cells_reverse(
    filtered.collect(),
    barcode_col=BARCODE,
    data_input=FILTERED_INPUT_DATA,
    top_cells=top_cells_df,
).with_columns(sample_name=pl.lit(sample_name))

# Fetch metrics

filtered_df_metrics_10000 = filtered_df_metrics.filter(
    pl.col("mean_filtered_reads_per_cell") == 10000
)
if not filtered_df_metrics_10000.is_empty():
    median_filtered_umis_per_cell_10000 = filtered_df_metrics_10000.select(
        "median_counts_per_cell"
    ).item()
    median_filtered_genes_per_cell_10000 = filtered_df_metrics_10000.select(
        "median_genes_per_cell"
    ).item()
else:
    median_filtered_umis_per_cell_10000 = None
    median_filtered_genes_per_cell_10000 = None
filtered_df_metrics_20000 = filtered_df_metrics.filter(
    pl.col("mean_filtered_reads_per_cell") == 20000
)
if not filtered_df_metrics_20000.is_empty():
    median_filtered_umis_per_cell_20000 = filtered_df_metrics_20000.select(
        "median_counts_per_cell"
    ).item()
    median_filtered_genes_per_cell_20000 = filtered_df_metrics_20000.select(
        "median_genes_per_cell"
    ).item()
else:
    median_filtered_umis_per_cell_20000 = None
    median_filtered_genes_per_cell_20000 = None


metrics.add("median_filtered_umis_per_cell_10000", median_filtered_umis_per_cell_10000)
metrics.add("median_filtered_umis_per_cell_20000", median_filtered_umis_per_cell_20000)
metrics.add(
    "median_filtered_genes_per_cell_10000", median_filtered_genes_per_cell_10000
)
metrics.add(
    "median_filtered_genes_per_cell_20000", median_filtered_genes_per_cell_20000
)
metrics.write_records("metrics.parquet")
del filtered
gc.collect()

filtered_df_metrics.write_csv(
    f"{sample_name}_{FILTERED_INPUT_DATA}_saturation_metrics.csv"
)

tenx_v31_filtered = REFERENCE_DATASETS["10x_v3.1_filtered"]

tenx_v4_filtered = REFERENCE_DATASETS["10x_v4_filtered"]


plot_comparative_saturation_curves(
    pl.concat([filtered_df_metrics, tenx_v31_filtered, tenx_v4_filtered]),
    sample_name=sample_name,
    input_type=FILTERED_INPUT_DATA,
)
