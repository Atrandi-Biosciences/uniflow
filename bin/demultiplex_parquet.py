#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

import polars as pl
from lib.common_const import (
    BARCODE,
    READ_NAME,
    SAMPLE_NAME,
)

library_id: str = sys.argv[1]
library_type: str = sys.argv[2]
barcode_mappings_path: Path = Path(sys.argv[3])
demultiplexing_sheet_path: Path = Path(sys.argv[4])
read1_path: Path = Path(sys.argv[5])
read2_path: Path = Path(sys.argv[6])
sample_name = sys.argv[7]

demultiplexing_sheet = pl.read_csv(demultiplexing_sheet_path)

# First check if the library needs to be demultiplexed or not

sample_df = (
    demultiplexing_sheet.filter(
        pl.col(SAMPLE_NAME) == sample_name,
        pl.col(BARCODE).is_not_null(),
    )
    .select(BARCODE)
    .unique()
)
if sample_df.select(BARCODE).is_empty():
    needs_demultiplexing = False
else:
    needs_demultiplexing = True
print(
    f"Library {library_id}/{library_type} needs demultiplexing: {needs_demultiplexing} for sample {sample_name}"
)
if needs_demultiplexing:
    barcode_mappings = pl.scan_parquet(barcode_mappings_path)
    subset = barcode_mappings.with_columns(
        barcode_a=pl.col(BARCODE).str.split("_").list.get(3)
    ).join(
        sample_df.lazy(),
        right_on=BARCODE,
        left_on="barcode_a",
        how="inner",
    )
    subset.sink_parquet(f"{sample_name}_barcode_mappings.parquet")
    subset.select(READ_NAME).unique().sink_csv(
        f"{sample_name}_reads.txt", include_header=False
    )
else:
    shutil.copy(
        barcode_mappings_path,
        f"{sample_name}_barcode_mappings.parquet",
        follow_symlinks=False,
    )
    empty_file = Path(f"{sample_name}_reads.txt")
    empty_file.touch()
