#!/usr/bin/env python

import sys
from pathlib import Path

import polars as pl
from lib.common_const import (
    LIBRARY_TYPE,
    MODALITY,
    PROCESS_NAME,
    SOURCE_ID,
)

experiment_id = sys.argv[1]
csv_metrics_list_path = sys.argv[2:]

metrics_list = []
for metrics_path in csv_metrics_list_path:
    metrics_path = Path(metrics_path)
    if metrics_path.suffix == ".csv":
        df = pl.scan_csv(
            metrics_path,
        )
    elif metrics_path.suffix == ".parquet":
        df = pl.scan_parquet(
            metrics_path,
        )
    else:
        continue

    metrics_list.append(df)

all_metrics_df = (
    pl.concat(metrics_list)
    .with_columns(pl.col("value").cast(pl.Float32))
    .sort([SOURCE_ID, LIBRARY_TYPE, MODALITY, PROCESS_NAME, "metric_key"])
)

all_metrics_df.sink_csv("metrics.csv")
