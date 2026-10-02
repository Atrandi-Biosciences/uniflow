#!/usr/bin/env python3
import sys
import polars as pl

source_id: str = sys.argv[1]
metrics_list: list = sys.argv[2:]

all_metrics = []
for metrics_path in metrics_list:
    all_metrics.append(pl.read_parquet(metrics_path))

merged_metrics = pl.concat(all_metrics)

merged_metrics.write_csv("metrics.csv")
