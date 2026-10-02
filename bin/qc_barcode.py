#!/usr/bin/env python3

import sys
from typing import Final

import polars as pl
from lib.common_const import (
    BARCODE,
    BARCODE_A,
    BARCODE_B,
    BARCODE_C,
    BARCODE_D,
    BARCODE_POSITION,
    HUMAN_READABLE_NAMES,
    LibraryType,
)
from lib.plotting.barcode import (
    plot_barcode_occurence_heatmap,
    plot_barcode_upset_plot,
    plot_per_base_substitution_rates,
    prepare_barcode_occurence_data,
)

# Declare inputs here
raw_bc_corrected_path = sys.argv[1]
substitution_parquet_path = sys.argv[2]
barcode_whitelist_path = sys.argv[3]
library_id = sys.argv[4]
library_type = sys.argv[5]

# Thresholds
top_fraction_used: Final[float] = 0.9
if library_type == HUMAN_READABLE_NAMES.get(LibraryType.RNA.value):
    WITH_UMI: bool = True
else:
    WITH_UMI: bool = False

ref_df = pl.read_csv(barcode_whitelist_path)
barcode_pool_size = int(ref_df.shape[0] / 4)
raw_bc_corrected = pl.scan_parquet(raw_bc_corrected_path).limit(1000000)
print("Upset Plots")

substitution_df = pl.read_parquet(substitution_parquet_path)
# edge case of error-free R2 (e.g., synthetic data)
if substitution_df.height:
    plot_per_base_substitution_rates(substitution_df, sample_name=library_id)


plot_barcode_upset_plot(
    raw_bc_corrected=raw_bc_corrected,
    ref_df=ref_df,
    library_id=library_id,
    library_type=library_type,
)

# Barcode pipetting heamap

filtered = (
    raw_bc_corrected.group_by(BARCODE_A, BARCODE_B, BARCODE_C, BARCODE_D)
    .agg(pl.len())
    .filter(
        pl.col(BARCODE_D).is_in(
            ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_D)[BARCODE]
        ),
        pl.col(BARCODE_C).is_in(
            ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_C)[BARCODE]
        ),
        pl.col(BARCODE_B).is_in(
            ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_B)[BARCODE]
        ),
        pl.col(BARCODE_A).is_in(
            ref_df.filter(pl.col(BARCODE_POSITION) == BARCODE_A)[BARCODE]
        ),
    )
    .sort("len", descending=True)
    .with_columns(pl.col("len").cum_sum().alias("cumsum"))
).collect()

total = filtered.select("cumsum").max().item()
cumsum_value = int(total * top_fraction_used)

plot_data = prepare_barcode_occurence_data(
    filtered.filter(pl.col("cumsum") < cumsum_value),
    reference_list=ref_df,
)

plot_barcode_occurence_heatmap(
    df=plot_data,
    cumsum_fraction=top_fraction_used,
    library_id=library_id,
    library_type=library_type,
    barcode_pool_size=barcode_pool_size,
)
