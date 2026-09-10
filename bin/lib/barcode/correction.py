#!/usr/bin/env python3

from typing import Tuple

import polars as pl
import polars_distance as pld
import numpy as np
from numpy.typing import NDArray

from lib.common_const import COUNTS, READS

# CONSTANTS


def get_per_base_substitution_rates(
    before_correction: pl.LazyFrame,
    after_correction: pl.LazyFrame,
    barcode_length: int,
    barcode_position: str,
    full_barcode_counts: pl.LazyFrame,
) -> pl.LazyFrame:
    splitted_before = (
        before_correction.join(
            full_barcode_counts, left_on="before", right_on=barcode_position
        )
        .with_columns(
            pl.col("before").str.split("").list.to_struct(upper_bound=barcode_length)
        )
        .unnest("before")
        .rename({COUNTS: "occurences"})
    )
    splitted_after = after_correction.with_columns(
        pl.col("after").str.split("").list.to_struct(upper_bound=barcode_length)
    ).unnest("after")
    omit_fields = [f"field_{i}" for i in range(barcode_length)]
    from_to = []
    total_reads = full_barcode_counts.select(pl.sum(COUNTS)).collect().item()
    for index, field in enumerate(omit_fields):
        joined = splitted_before.join(
            splitted_after,
            on=splitted_before.select(pl.exclude(field, "occurences")).columns,
            how="inner",
        )
        substitutions = (
            joined.select(field, f"{field}_right", "occurences")
            .group_by(pl.all())
            .agg(pl.len().alias("per_substitution_per_index_count"))
            .with_columns(pl.lit(index).alias("barcode_index"))
            .rename({field: "mutated_base", f"{field}_right": "original_base"})
            .filter(pl.col("original_base") != pl.col("mutated_base"))
            .with_columns(
                total_counts=pl.col("per_substitution_per_index_count")
                * pl.col("occurences"),
            )
        )
        from_to.append(substitutions)
    substitutions_df = (
        pl.concat(from_to)
        .with_columns(
            total_reads=pl.lit(total_reads),
        )
        .with_columns(
            (pl.col("total_counts") / (pl.col("total_reads") * barcode_length)).alias(
                "substitution_frequency"
            ),
        )
        .with_columns(
            barcode_position=pl.lit(barcode_position),
            feature_id=pl.col("original_base") + ">" + pl.col("mutated_base"),
        )
    )
    return substitutions_df


def correct_barcodes_pl(
    barcodes_df: pl.DataFrame,
    barcode_subset_df: pl.DataFrame,
    hamming_distance: int,
    barcode_col_name: str,
    whitelist_col_name: str,
) -> Tuple[dict, pl.DataFrame]:
    """Corrects barcodes using a subset based on join_asof from polars.
    Uses both forward and backward strategy to dinf the closest barcode

    Args:
        barcodes_df (pl.DataFrame): All barcodes with their respective counts
        barcode_subset_df (pl.DataFrame): Barcode reference used to correct
        hamming_distance (int): Max hamming distance allowed
        mapped_barcodes (dict): Dict of mapped barcodes
        pl.DataFrame: Corrected barcodes including the original and whitelist columns

    Returns:
        tuple[pl.DataFrame, int]: The corrected version of the input barcodes_df, number of corrected barcodes
    """
    print("Correcting barcodes")
    corrected_barcodes_pl = pl.DataFrame(
        schema={
            barcode_col_name: pl.String,
            COUNTS: pl.UInt32,
            whitelist_col_name: pl.String,
        }
    )
    current_barcodes = pl.DataFrame(
        schema={
            barcode_col_name: pl.String,
            COUNTS: pl.UInt32,
            whitelist_col_name: pl.String,
        }
    )
    current_barcodes_to_correct = barcodes_df.shape[0]
    last_iteration_barcodes_to_correct = 0
    unknown_barcodes = barcodes_df.filter(
        ~pl.col(barcode_col_name).is_in(barcode_subset_df[whitelist_col_name])
    )
    n_iterations = 0
    while (
        current_barcodes_to_correct > 0
        and current_barcodes_to_correct != last_iteration_barcodes_to_correct
    ):
        methods = ["backward", "forward"]
        for method in methods:
            current_barcodes = (
                (
                    unknown_barcodes.filter(
                        (
                            ~pl.col(barcode_col_name).is_in(
                                corrected_barcodes_pl[barcode_col_name]
                            )
                        )
                    )
                    .sort(barcode_col_name)
                    .join_asof(
                        barcode_subset_df.sort(whitelist_col_name),
                        left_on=barcode_col_name,
                        right_on=whitelist_col_name,
                        strategy=method,  # type: ignore
                    )
                )
                .filter(~pl.col(whitelist_col_name).is_null())
                .with_columns(
                    pld.col(barcode_col_name)
                    .dist_str.hamming(pl.col(whitelist_col_name))
                    .cast(pl.UInt32)
                    .alias("hamming_distance")
                )
                .filter(pl.col("hamming_distance") <= hamming_distance)
                .drop("hamming_distance")
            )
            corrected_barcodes_pl = pl.concat(
                [
                    corrected_barcodes_pl,
                    current_barcodes,
                ]
            )
        current_barcodes_to_correct = current_barcodes.shape[0]
        barcode_subset_df = barcode_subset_df.filter(
            ~pl.col(whitelist_col_name).is_in(corrected_barcodes_pl[whitelist_col_name])
        )
        n_iterations += 1
    print(f"Corrected barcodes in {n_iterations} iterations")
    print(f"Number of uncorrected barcodes: {current_barcodes.shape[0]}")
    mapped_barcodes = dict(
        corrected_barcodes_pl.select(barcode_col_name, whitelist_col_name).iter_rows()  # type: ignore
    )
    print("Barcodes corrected")

    return mapped_barcodes, corrected_barcodes_pl.select(
        barcode_col_name, whitelist_col_name
    )


def get_motif_freq(one_column_df: pl.DataFrame, size: int) -> NDArray | None:
    """Get motif frequency of DNA bases

    Args:
        one_column_df (pl.DataFrame): On column df of the sequence we want the motif from
        size (int): length of the sequence

    Returns:
        NDArray | None: array struct for plotting
    """
    A = []
    C = []
    G = []
    T = []
    name = one_column_df.columns[0]
    splitted = (
        one_column_df.filter(
            pl.col(name).is_not_null(),
            pl.col(name).str.len_chars() == 8,
            ~pl.col(name).str.contains("N"),
        )
        .with_columns(pl.col(name).str.split_exact(by="", n=(size - 1)))
        .unnest(name)
    )
    for name in splitted.columns:
        aggregated = splitted.select(name).group_by(name).agg(pl.len().alias(READS))
        A.append(aggregated.filter(pl.col(name) == "A")[READS].item())
        C.append(aggregated.filter(pl.col(name) == "C")[READS].item())
        G.append(aggregated.filter(pl.col(name) == "G")[READS].item())
        T.append(aggregated.filter(pl.col(name) == "T")[READS].item())
    freq_list = np.array([np.array(A), np.array(C), np.array(G), np.array(T)])
    return freq_list
