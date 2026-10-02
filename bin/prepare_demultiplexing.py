#!/usr/bin/env python3

import sys
from pathlib import Path

import polars as pl
from lib.common_const import (
    BARCODE,
    BARCODE_POSITION,
    DEMULTIPLEXING_INDICES,
    FASTQ_1,
    FASTQ_2,
    INDEX,
    LIBRARY_ID,
    LIBRARY_TYPE,
    SAMPLE_NAME,
)

samplesheet_csv_path: Path = Path(sys.argv[1])
barcode_whitelist_path: Path = Path(sys.argv[2])
# TODO: Implement user input, currently hardcoded to A
barcode_position_for_demultiplexing: str = sys.argv[3]

sample_sheet = pl.read_csv(samplesheet_csv_path)

barcode_plate = pl.read_csv(barcode_whitelist_path).filter(
    pl.col(BARCODE_POSITION) == barcode_position_for_demultiplexing
)

demultiplexing_sheet = sample_sheet.with_columns(
    library_id=pl.col(FASTQ_2).str.extract(
        "^.*/(.*)_R2*"
    )  # take note of the "." at the end
)  # .select(  # TODO: This is very brittle, we should have a more robust way to handle it
#     SAMPLE_NAME, demultiplexing_indices, LIBRARY_ID, LIBRARY_TYPE
# )


if DEMULTIPLEXING_INDICES in demultiplexing_sheet.columns:
    fully_melted = demultiplexing_sheet.with_columns(
        pl.col(DEMULTIPLEXING_INDICES).str.split(by=";")
    ).explode(DEMULTIPLEXING_INDICES)
    sample_mapping = barcode_plate.join(
        fully_melted, right_on=DEMULTIPLEXING_INDICES, left_on=INDEX, how="right"
    )
else:
    sample_mapping = (
        demultiplexing_sheet.with_columns(
            pl.lit(None).alias(BARCODE_POSITION),
            pl.lit(None).alias(BARCODE),
            pl.lit(None).alias(DEMULTIPLEXING_INDICES),
        )
    ).select(
        BARCODE_POSITION,
        BARCODE,
        SAMPLE_NAME,
        FASTQ_1,
        FASTQ_2,
        LIBRARY_TYPE,
        LIBRARY_ID,
        DEMULTIPLEXING_INDICES,
    )


sample_mapping.write_csv("sample_demultiplexing_sheet.csv")

sample_sheet.write_csv("sample_sheet.csv")
