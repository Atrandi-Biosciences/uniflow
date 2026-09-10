#!/usr/bin/env python3

import sys
import gc
import polars as pl
import os
from pathlib import Path
from typing import Final


from lib.barcode.correction import correct_barcodes_pl, get_per_base_substitution_rates
from lib.dataframe_toolkit import get_most_freq_df
from lib.metrics.metrics import MetricsRecord
from lib.common_const import (
    BARCODE_D,
    BARCODE_C,
    BARCODE_B,
    BARCODE_A,
    BARCODE,
    READ_NAME,
    COUNTS,
    PRIMER,
    BARCODE_POSITION,
    INDEX,
    Modality,
    LibraryType,
)

jemalloc_conf = "dirty_decay_ms:500,muzzy_decay_ms:-1"
if os.environ.get("POLARS_THP") == "1":
    jemalloc_conf += ",thp:always,metadata_thp:always"
if override := os.environ.get("_RJEM_MALLOC_CONF"):
    jemalloc_conf += "," + override
os.environ["_RJEM_MALLOC_CONF"] = jemalloc_conf

# Declare inputs here
barcode_parquet_path: Path = Path(sys.argv[1])
barcode_whitelist_path: Path = Path(sys.argv[2])
library_type: str = sys.argv[3]
sample_name: str = sys.argv[4]
HAMMING_DISTANCE: Final[int] = 1

metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType(library_type),
    modality=Modality.BARCODE,
)

ref_df = pl.read_csv(barcode_whitelist_path)
barcodes: pl.LazyFrame = pl.scan_parquet(barcode_parquet_path)

df_indexes = list(barcodes.drop(READ_NAME).collect_schema().keys())

if PRIMER in df_indexes:
    additional_ref = (
        get_most_freq_df(barcodes, key=PRIMER)
        .rename({PRIMER: BARCODE})
        .collect()
        .with_columns(pl.lit(PRIMER).alias(BARCODE_POSITION), pl.lit(None).alias(INDEX))
        .select(BARCODE_POSITION, BARCODE, INDEX)
    )
    # Change this once we have a better way of doing it.
    additional_ref = pl.DataFrame(
        {BARCODE_POSITION: [PRIMER], BARCODE: ["GACTTGAGTGGCTGTCGG"], INDEX: [None]}
    )
    ref_df = pl.concat([ref_df, additional_ref])

corrected_mapping = {i: {} for i in df_indexes}

raw_reads = barcodes.collect(engine="streaming")
substitutions_list = []
for position_to_correct in corrected_mapping:
    # Barcode correction
    per_position_barcode_counts = (
        raw_reads.select(position_to_correct)
        .group_by(position_to_correct)
        .agg(pl.len().alias(COUNTS))
    )
    mapped, mapped_df = correct_barcodes_pl(
        per_position_barcode_counts,
        ref_df.filter(pl.col(BARCODE_POSITION) == position_to_correct).select(BARCODE),
        hamming_distance=HAMMING_DISTANCE,
        barcode_col_name=position_to_correct,
        whitelist_col_name=BARCODE,
    )
    corrected_mapping[position_to_correct] = mapped

    if position_to_correct in [BARCODE_A, BARCODE_B, BARCODE_C, BARCODE_D]:
        substitutions = get_per_base_substitution_rates(
            before_correction=mapped_df.lazy()
            .select(position_to_correct)
            .rename({position_to_correct: "before"}),
            after_correction=mapped_df.lazy()
            .select(BARCODE)
            .rename({BARCODE: "after"}),
            barcode_length=8,
            barcode_position=position_to_correct,
            full_barcode_counts=per_position_barcode_counts.lazy(),
        )
        substitutions_list.append(substitutions)


all_substitutions = pl.concat(substitutions_list)
all_substitutions.sink_parquet(f"{sample_name}_substitutions.parquet")

total_reads = raw_reads.shape[0]
metrics.add("total_reads", total_reads)

corrected_columns = list(corrected_mapping)

expressions = [
    pl.col(col_name).replace_strict(
        corrected_mapping[col_name], default=pl.col(col_name)
    )
    for col_name in corrected_columns
]

raw_bc_corrected = raw_reads.with_columns(*expressions).lazy()


del raw_reads
gc.collect()

raw_bc_corrected.sink_parquet(f"{sample_name}_corrected_barcode_reads.parquet")


barcode_mapping = raw_bc_corrected.filter(
    pl.col(BARCODE_D).is_in(
        ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_D)[BARCODE]
    )
    & pl.col(BARCODE_C).is_in(
        ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_C)[BARCODE]
    )
    & pl.col(BARCODE_B).is_in(
        ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_B)[BARCODE]
    )
    & pl.col(BARCODE_A).is_in(
        ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_A)[BARCODE]
    )
)

del raw_bc_corrected
gc.collect()

# We want to keep the UIM column when running RNA libraries
full_barcode_mapping = barcode_mapping.with_columns(
    barcode=pl.concat_str(BARCODE_D, BARCODE_C, BARCODE_B, BARCODE_A, separator="_")
).drop(
    BARCODE_A, BARCODE_B, BARCODE_C, BARCODE_D, strict=False
)  # Might be dangerous, find a better way

# This is equivalent to df.shape[0] for a lazy frame
corrected_filtered_reads = barcode_mapping.select(pl.len()).collect().item()

metrics.add("total_reads_with_valid_full_barcode", corrected_filtered_reads)
metrics.add(
    "total_reads_without_valid_full_barcode", total_reads - corrected_filtered_reads
)
fraction_reads_with_valid_full_barcode = corrected_filtered_reads / total_reads
metrics.add(
    "fraction_reads_with_valid_full_barcode",
    fraction_reads_with_valid_full_barcode,
)
metrics.add(
    "fraction_reads_without_valid_full_barcode",
    1 - fraction_reads_with_valid_full_barcode,
)
full_barcode_mapping.sink_parquet(f"{sample_name}_barcodes_mapping.parquet")

# Write out metrics

metrics.write_records("metrics.parquet")
