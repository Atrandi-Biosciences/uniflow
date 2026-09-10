#!/usr/bin/env python3
import sys
import polars as pl

source_id: str = sys.argv[1]
status_list: list = sys.argv[2:]

all_status = []

for status_path in status_list:
    all_status.append(pl.scan_parquet(status_path))

all_status_df = pl.concat(all_status)
all_status_df.sink_csv("status.csv")
