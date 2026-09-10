#!/usr/bin/env python3

import polars as pl

import sys

sample_name = sys.argv[1]
library_type = sys.argv[2]
input_parquets = sys.argv[3].split(" ")

pl.scan_parquet(input_parquets).sink_parquet(
    f"{sample_name}_{library_type}_barcode_mappings.parquet"
)
