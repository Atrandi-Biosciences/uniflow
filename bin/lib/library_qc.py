import polars as pl
import numpy as np
import plotnine as gg
from lib.common_const import MAPPING_QUALITY, COUNTS
from lib.bam_tools.utils import (
    GENE_SYMBOL_TAG,
    CORRECTED_UMI_TAG,
    STAR_BAD_UMI,
    STAR_NO_GENE,
    STAR_MIN_MAPPING_QUALITY,
)
from typing import Final

# constants
N_GENES: Final[str] = "n_genes"
MEDIAN_COUNTS_PER_CELL: Final[str] = "median_counts_per_cell"
MEDIAN_GENES_PER_CELL: Final[str] = "median_genes_per_cell"
SATURATION: Final[str] = "saturation"


# Use this for full lazyframe implementation
def subsample_lf(
    lf: pl.LazyFrame, total_reads: int, n_subsampled_reads: int
) -> pl.LazyFrame:
    """Subsamples a LazyFrame to n_subsampled_reads by randomly selecting indices."""
    indices = sorted(
        np.random.choice(total_reads, size=n_subsampled_reads, replace=False)
    )
    subsampled_reads_lf = lf.select(pl.col("*").gather(indices))
    return subsampled_reads_lf


def generate_round_descending_range(max_value: int, step: int) -> list[int]:
    """Generates a descending list starting from max_value, followed by the largest multiples
    of step that are <= max_value, stepping down by step until 0.

    The list always includes max_value (if > 0) and ends with 0. If max_value is already a
    multiple of step, it is included only once.

    Args:
        max_value: The maximum starting value (non-negative integer).
        step: The positive step size for multiples (positive integer).

    Returns:
        A list of integers in descending order.

    Raises:
        ValueError: If max_value < 0 or step <= 0.
    """
    if max_value < 0:
        raise ValueError("max_value must be non-negative")
    if step <= 0:
        raise ValueError("step must be positive")

    if max_value == 0:
        return [0]

    largest_multiple = (max_value // step) * step
    multiples = list(range(largest_multiple, step, -step))

    return multiples


def get_saturation_metrics_df(df_raw: pl.DataFrame, barcode_col: str, n_reads: int):
    """Generate and retrieve saturation metrics from a polars dataframe at a specific read depth.

    Args:
        df_raw (pl.DataFrame): raw polars DataFrame containing read data
        barcode_col (str): name of the column containing barcode information
        n_reads (int): total number of reads in the dataset

    Returns:
        tuple: containing median counts per cell, median genes per cell, and saturation rate
    """
    raw_aggregated = df_raw.group_by(
        barcode_col, GENE_SYMBOL_TAG, CORRECTED_UMI_TAG, MAPPING_QUALITY
    ).agg(pl.len().alias(COUNTS))

    df_aggregated = raw_aggregated.filter(
        pl.col(MAPPING_QUALITY) == STAR_MIN_MAPPING_QUALITY
    ).filter(
        pl.col(GENE_SYMBOL_TAG) != STAR_NO_GENE,
        pl.col(GENE_SYMBOL_TAG).is_not_null(),
        pl.col(CORRECTED_UMI_TAG) != STAR_BAD_UMI,
    )

    median_counts_per_cell = (
        df_aggregated.group_by(barcode_col, GENE_SYMBOL_TAG)
        .agg(pl.len().alias(COUNTS))
        .group_by(barcode_col)
        .agg(pl.sum(COUNTS))
        .select(pl.median(COUNTS))
        .item()
    )
    median_genes_per_cell = (
        df_aggregated.group_by(barcode_col, GENE_SYMBOL_TAG)
        .agg(pl.len().alias(COUNTS))
        .group_by(barcode_col)
        .agg(pl.len().alias(N_GENES))
        .select(pl.median(N_GENES))
        .item()
    )
    saturation_rate = 1 - (raw_aggregated.shape[0] / n_reads)
    return (median_counts_per_cell, median_genes_per_cell, saturation_rate)


def get_saturation_curve_cells_reverse(
    full_read: pl.DataFrame,
    data_input: str,
    barcode_col: str,
    top_cells: pl.DataFrame,
    step: int = 5000,
) -> pl.DataFrame:
    """Generate saturation curves by subsampling reads per cell in descending order.

    Args:
        full_read (pl.DataFrame): DataFrame containing full read data
        data_input (str): Type of input data (e.g., "raw" or "filtered")
        barcode_col (str): Name of the column containing barcode information
        top_cells (pl.DataFrame): DataFrame containing top cells information
        step (int, optional): Step size for subsampling reads per cell. Defaults to 5000.

    Returns:
        pl.DataFrame: DataFrame containing saturation metrics
    """
    n_cells = top_cells.shape[0]
    total_reads = full_read.shape[0]
    reads_per_cell = generate_round_descending_range(total_reads // n_cells, step)
    reads_per_cell.extend([2500, 1000, 500])
    saturation_metrics = {
        MEDIAN_COUNTS_PER_CELL: [0],
        MEDIAN_GENES_PER_CELL: [0],
        f"mean_{data_input}_reads_per_cell": [0],
        SATURATION: [0.0],
    }
    median_counts_per_cell, median_genes_per_cell, saturation_rate = (
        get_saturation_metrics_df(
            full_read.join(top_cells, on=barcode_col),
            barcode_col=barcode_col,
            n_reads=total_reads,
        )
    )
    saturation_metrics[f"mean_{data_input}_reads_per_cell"].append(
        int(total_reads / n_cells)
    )
    saturation_metrics[MEDIAN_COUNTS_PER_CELL].append(median_counts_per_cell)
    saturation_metrics[MEDIAN_GENES_PER_CELL].append(median_genes_per_cell)
    saturation_metrics[SATURATION].append(saturation_rate)
    for i in reads_per_cell:
        subsampled_reads = n_cells * i
        if subsampled_reads > full_read.shape[0]:
            print(
                f"Skipping subsample size {subsampled_reads} as it exceeds full read count {full_read.shape[0]}"
            )
            continue
        full_read = full_read.sample(subsampled_reads, seed=42)
        median_counts_per_cell, median_genes_per_cell, sat_rate = (
            get_saturation_metrics_df(
                full_read.join(top_cells, on=barcode_col),
                barcode_col=barcode_col,
                n_reads=subsampled_reads,
            )
        )
        saturation_metrics[f"mean_{data_input}_reads_per_cell"].append(i)
        saturation_metrics[MEDIAN_COUNTS_PER_CELL].append(median_counts_per_cell)
        saturation_metrics[MEDIAN_GENES_PER_CELL].append(median_genes_per_cell)
        saturation_metrics[SATURATION].append(sat_rate)

    return df_from_saturation_dict(input_dict=saturation_metrics, data_input=data_input)


def df_from_saturation_dict(input_dict: dict, data_input: str) -> pl.DataFrame:
    """Convert from a dict metrics to a df

    Args:
        input_dict (dict): a simple key,value pairs of saturation metrics
        data_input (str): raw or filtered

    Returns:
        pl.DataFrame: df version of the same data
    """
    df = pl.DataFrame(
        input_dict,
        schema={
            f"mean_{data_input}_reads_per_cell": pl.UInt32,
            MEDIAN_COUNTS_PER_CELL: pl.Float32,
            MEDIAN_GENES_PER_CELL: pl.Float32,
            SATURATION: pl.Float32,
        },
    )
    return df


def plot_saturation_curves(
    metrics_df: pl.DataFrame, sample_name: str, input_type: str
) -> None:
    """Plot 3 saturation curves

    Args:
        metrics_df (pl.DataFrame): metrics for saturation curves
        sample_name (str): name of the sample
        input_type (str): raw or filtered
    """
    plot = (
        gg.ggplot(metrics_df.melt(id_vars=f"mean_{input_type}_reads_per_cell"))
        + gg.aes(x=f"mean_{input_type}_reads_per_cell", y="value")
        + gg.facet_wrap("variable", scales="free")
        + gg.geom_line()
        + gg.theme(figure_size=(10, 4))
        + gg.ggtitle(f"Saturation curves for {sample_name}", subtitle=input_type)
    )
    plot.save(f"{sample_name}_{input_type}_saturation.png")


def plot_comparative_saturation_curves(
    metrics_df: pl.DataFrame, sample_name: str, input_type: str
) -> None:
    """Plot 3 saturation curves

    Args:
        metrics_df (pl.DataFrame): metrics for saturation curves
        sample_name (str): name of the sample
        input_type (str): raw or filtered
    """
    plot = (
        gg.ggplot(
            metrics_df.melt(
                id_vars=[f"mean_{input_type}_reads_per_cell", "sample_name"]
            )
        )
        + gg.aes(x=f"mean_{input_type}_reads_per_cell", y="value", color="sample_name")
        + gg.facet_wrap("variable", scales="free")
        + gg.geom_line()
        + gg.theme(figure_size=(10, 4))
        + gg.ggtitle(f"Saturation curves for {sample_name}", subtitle=input_type)
    )
    plot.save(f"{input_type}_saturation_mqc.png")


REFERENCE_DATASETS = {
    "ref4": pl.DataFrame(
        {
            "mean_filtered_reads_per_cell": [
                0,
                500,
                1000,
                2500,
                10000,
                15000,
                20000,
                22403,
            ],
            "median_counts_per_cell": [
                0.0,
                453.0,
                886.0,
                2104.0,
                6611.0,
                8662.5,
                10238.0,
                10868.5,
            ],
            "median_genes_per_cell": [
                0.0,
                335.0,
                579.0,
                1109.5,
                2340.0,
                2715.0,
                2963.0,
                3056.0,
            ],
            "saturation": [
                0.0,
                0.0182736,
                0.0361918,
                0.0862732,
                0.27820757,
                0.36879548,
                0.4400414,
                0.46910378,
            ],
        },
        schema={
            "mean_filtered_reads_per_cell": pl.UInt32,
            "median_counts_per_cell": pl.Float32,
            "median_genes_per_cell": pl.Float32,
            "saturation": pl.Float32,
        },
    ).with_columns(sample_name=pl.lit("ref4")),
    "ref3": pl.DataFrame(
        {
            "mean_filtered_reads_per_cell": [
                0,
                500,
                1000,
                2500,
                10000,
                15000,
                20000,
                25000,
                30000,
                35000,
                40000,
                45000,
                47899,
            ],
            "median_counts_per_cell": [
                0.0,
                461.0,
                897.0,
                2094.0,
                6110.0,
                7698.0,
                8789.0,
                9600.0,
                10203.0,
                10664.0,
                11035.0,
                11331.0,
                11475.0,
            ],
            "median_genes_per_cell": [
                0.0,
                350.0,
                607.0,
                1174.0,
                2483.0,
                2863.0,
                3101.0,
                3258.0,
                3372.0,
                3457.0,
                3527.0,
                3578.0,
                3605.0,
            ],
            "saturation": [
                0.0,
                0.024662351,
                0.048454188,
                0.11407432,
                0.35179967,
                0.45536134,
                0.532512,
                0.5915316,
                0.6379222,
                0.6750953,
                0.7054708,
                0.7306858,
                0.7434344,
            ],
        },
        schema={
            "mean_filtered_reads_per_cell": pl.UInt32,
            "median_counts_per_cell": pl.Float32,
            "median_genes_per_cell": pl.Float32,
            "saturation": pl.Float32,
        },
    ).with_columns(sample_name=pl.lit("ref3")),
}
