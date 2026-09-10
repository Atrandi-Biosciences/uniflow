#!/usr/bin/env python3

import polars as pl
import sys

from lib.common_const import (
    BARCODE_D,
    BARCODE_C,
    BARCODE_B,
    BARCODE_A,
    PRIMER,
    UMI,
    READ_NAME,
)

# Declare inputs here
barcode_fastq_csv_path = sys.argv[1]
library_type = sys.argv[2]
sample_name = sys.argv[3]

if library_type == "RNA":
    WITH_UMI: bool = True
else:
    WITH_UMI: bool = False
HAMMING_DISTANCE = 1
barcode_length: int = 70 if WITH_UMI else 48


if WITH_UMI:
    barcodes: pl.LazyFrame = pl.scan_csv(
        barcode_fastq_csv_path,
        schema={
            READ_NAME: pl.String,
            BARCODE_D: pl.String,
            BARCODE_C: pl.String,
            BARCODE_B: pl.String,
            BARCODE_A: pl.String,
            PRIMER: pl.String,
            UMI: pl.String,
        },
        has_header=False,
    )
else:
    barcodes: pl.LazyFrame = pl.scan_csv(
        barcode_fastq_csv_path,
        schema={
            READ_NAME: pl.String,
            BARCODE_D: pl.String,
            BARCODE_C: pl.String,
            BARCODE_B: pl.String,
            BARCODE_A: pl.String,
        },
        has_header=False,
    )

barcodes.sink_parquet(f"{sample_name}_barcodes.parquet")
