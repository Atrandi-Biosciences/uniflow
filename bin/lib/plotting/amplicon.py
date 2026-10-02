import plotnine as gg
import polars as pl
import math
import numpy as np
from lib.common_const import CHROM


def square_subplot_shape(n_features: int) -> tuple[int, int]:
    """
    Return (n_rows, n_cols) that is as square as possible and large
    enough to hold n_features subplots.
    """
    if n_features <= 0:
        return (1, 1)

    # Start with the integer closest to the square root
    cols = math.ceil(math.sqrt(n_features))
    rows = math.ceil(n_features / cols)

    # Sometimes trimming one column makes it more square
    if cols > 1 and (cols - 1) * rows >= n_features:
        cols -= 1

    return rows, cols


def plot_amplicon_anchors(
    positions: pl.DataFrame, sample_name: str, title: str, prefix: str
):
    n_features = positions.select(CHROM).unique().shape[0]
    n_row, n_col = square_subplot_shape(n_features)
    plot = (
        gg.ggplot(positions)
        + gg.aes(x="start", y="reads", color="mean_mapq")
        + gg.geom_point(size=0.5)
        + gg.facet_wrap(CHROM, ncol=n_col, scales="free_x")
        + gg.theme(figure_size=(1.5 + 2 * n_col, 1 + 1.5 * n_row))
        + gg.scale_color_continuous("turbo")
        + gg.ggtitle(f"{title}\n{sample_name}")
        + gg.scale_y_log10()
        + gg.xlab("Alignment start position")
    )
    plot.save(f"{prefix}_alignment_distribution.png")


def get_brackets_from_breaks(breaks):
    """Generate brackets for a list of breaks. Add 0-first and last-

    Args:
        breaks (list[int]): List of break points.

    Returns:
        list[str]: List of bracket labels.
    """
    brackets = []
    brackets.append(f"0-{breaks[0]}")
    for i in range(len(breaks) - 1):
        brackets.append(f"{breaks[i]}-{breaks[i + 1]}")
    brackets.append(f"{breaks[-1]}-")
    return brackets


def complete_amplicon_binned_df(binned: pl.DataFrame) -> pl.DataFrame:
    """Generate missing full cross join of binned amplicon data counts

    Args:
        binned (pl.DataFrame): Binned amplicon data.

    Returns:
        pl.DataFrame: Complete binned amplicon data with missing combinations filled.
    """
    plot_data = (
        binned.group_by(CHROM, "breaks")
        .agg(pl.len().alias("n_cells"))
        .with_columns(pl.col("breaks"))
    )
    full_data = (
        binned.select(CHROM)
        .unique()
        .join(binned.select("breaks").unique(), how="cross")
    )
    full_plot_data = plot_data.join(
        full_data, on=[CHROM, "breaks"], how="right"
    ).fill_null(0)
    return full_plot_data


def plot_amplicon_performance(
    filtered_counts: pl.DataFrame,
    breaks: list[int],
    sample_name: str,
    min_threshold: int,
):
    """Plot amplicon performance.

    Args:
        filtered_counts (pl.DataFrame): Filtered amplicon counts. Required columns: "chrom", "barcode", "reads".
        breaks (list[int]): List of break points.
        sample_name (str): Sample name.
        min_threshold (int): Minimum threshold for reads.
    """
    min_bin_over_threshold = np.argmax(np.array(breaks) >= min_threshold) + 1.5
    brackets = get_brackets_from_breaks(breaks)
    binned = filtered_counts.with_columns(
        pl.col("reads").cut(breaks=breaks, labels=brackets).alias("breaks")
    )
    full_plot_data = complete_amplicon_binned_df(binned)
    plot = (
        gg.ggplot(full_plot_data, gg.aes(y="breaks", x=CHROM, fill="n_cells"))
        + gg.geom_tile()
        + gg.coord_fixed()
        + gg.theme(axis_text_x=gg.element_text(rotation=45, ha="right"))
        + gg.scale_y_discrete(limits=brackets)
        + gg.scale_fill_continuous(na_value="#440154")
        + gg.labs(
            x="Amplicon",
            y="Reads",
            title=f"{sample_name}: Number of cells per bin",
            subtitle=f"Total cells: {filtered_counts.select('barcode').unique().shape[0]}",
        )
        + gg.geom_hline(
            yintercept=min_bin_over_threshold, linetype="dashed", color="white", size=1
        )
    )
    plot.save("amplicon_performance.png")
