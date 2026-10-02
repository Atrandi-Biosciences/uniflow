#!/usr/bin/env python3

import polars as pl
from lib.plotting.amplicon import plot_amplicon_anchors, plot_amplicon_performance
from lib.plotting.common import (
    plot_generic_upset_plot,
    export_upset_data_csv,
)
from lib.common_const import FRACTION, READS, Modality, label_for, LibraryType
from lib.metrics.metrics import MetricsRecord

import sys

sample_name = sys.argv[1]
position_parquet = sys.argv[2]
upset_data_parquet = sys.argv[3]
filtered_counts_parquet = sys.argv[4]
top_cells_parquet = sys.argv[5]

positions = pl.read_parquet(position_parquet)
bool_counts_df = pl.read_parquet(upset_data_parquet)
counts = pl.read_parquet(filtered_counts_parquet)
top_cells = pl.read_parquet(top_cells_parquet)
filtered_counts = counts.join(top_cells, on="barcode", how="inner")

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.AMPLICON,
)
plot_amplicon_anchors(
    positions=positions,
    sample_name=sample_name,
    title="Amplicon mapping distributions",
    prefix="amplicon",
)
plot_amplicon_performance(
    filtered_counts=filtered_counts,
    breaks=[0, 1, 10, 100, 1000, 2000, 5000, 10000],
    sample_name=sample_name,
    min_threshold=80,
)
plot_generic_upset_plot(
    bool_count=bool_counts_df,
    sample_name=sample_name,
    title=f"{label_for(Modality.AMPLICON.value)} construct integrity",
    modularity=Modality.AMPLICON.value,
)

amplicon_upset_csv = export_upset_data_csv(
    bool_lf=pl.scan_parquet(upset_data_parquet),
    sample_name=sample_name,
    prefix="amplicon",
)

amplicon_upset_metrics = pl.read_csv(amplicon_upset_csv)

usable_reads = amplicon_upset_metrics.filter(
    pl.col("has_barcode")
    & pl.col("is_anchored")
    & pl.col("mapped")
    & pl.col("high_quality_mapped")
)

fraction_usable_reads = usable_reads.select(FRACTION).item()
total_usable_reads = usable_reads.select(READS).item()
metrics.add("fraction_usable_reads", fraction_usable_reads)
metrics.add("total_usable_reads", total_usable_reads)
metrics.write_records("metrics.parquet")
