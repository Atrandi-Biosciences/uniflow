#!/usr/bin/env python3

import sys
import polars as pl

from lib.common_const import READ_NAME

bam_parquet_path = sys.argv[1]
barcode_mappings_path = sys.argv[2]
sample_name = sys.argv[3]
library_type = sys.argv[4]


bam_lf = pl.scan_parquet(bam_parquet_path)
barcode_mappings_lf = pl.scan_parquet(barcode_mappings_path)

merged = bam_lf.join(
    barcode_mappings_lf,
    on=READ_NAME,
    how="left",
)
merged.sink_parquet(f"full_read_{library_type}.parquet")
