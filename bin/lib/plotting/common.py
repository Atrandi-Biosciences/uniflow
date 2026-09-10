from typing import List

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from lib.common_const import FRACTION, READS, label_for
from matplotlib import rcParams
from upsetplot import UpSet, from_memberships


def export_upset_data_csv(
    bool_lf: pl.LazyFrame, sample_name: str, prefix: str, min_fraction: float = 0.01
) -> str:
    out_path = f"{prefix}_upset_data.csv"
    bool_lf.with_columns(
        (pl.col(READS) / pl.sum(READS)).alias(FRACTION), sample_name=pl.lit(sample_name)
    ).sort(FRACTION, descending=True).sink_csv(out_path)
    return out_path


def convert_df_to_upset_data(df: pl.DataFrame) -> tuple[list[list], list]:
    """Convert a boolean dataframe to an upsetplot compatible format

    Args:
        df (pl.DataFrame): boolean df of whitelisted barcodes

    Returns:
        Tuple[list, list]: List of groups, list of read counts
    """
    col_indexes = df.columns
    col_indexes.remove(READS)
    final_list: list[list] = []
    read_counts = []
    for i in df.iter_rows():
        value = i[-1]
        current_list = []
        for index, j in enumerate(i):
            if index > len(col_indexes):
                if i[index] > 0:
                    continue
            if j is True:
                current_list.append(col_indexes[index])

        final_list.append(current_list)
        read_counts.append(value)
    return final_list, read_counts


def convert_membership_names(memberships: list[list]) -> list[list]:
    """Tries to translate variables to human readable names for upset memberships
    Defaults to the variable name if HUMAN_READABLE_NAMES does not have the required entry

    Args:
        memberships (list[list]): variable memberships

    Returns:
        list[list]: Human readable memberships
    """
    new_memberships = []
    for group in memberships:
        new_group = []
        for element in group:
            new_element = label_for(key=element)
            new_group.append(new_element)
        new_memberships.append(new_group)
    return new_memberships


def plot_percell_variant_abundance(
    sites: pl.DataFrame,
    sample_name: str,
    output_path: str,
    ncallers_high: int = 2,
    multicell_min: int = 2,
    variant_class: str = "",
) -> None:
    """Per-site per-cell-vs-abundance scatter

    Each point is one called site; x is caller agreement (`n_callers`, null
    treated as 0) and y is cell support (`cells_supporting_alt`). The dashed
    guides at `ncallers_high` and `multicell_min` split the plot into quadrants;
    the top-right corner ; called by enough callers AND backed by enough cells ;
    is the high-confidence set.

    Args:
        sites: one row per (feature, pos_local, ref, alt) with `n_callers` and
            `cells_supporting_alt` columns.
        sample_name: sample label for the title.
        output_path: png path to write.
        ncallers_high: caller-agreement boundary (matches consensus high-confidence).
        multicell_min: cell-support boundary (matches multicell_support).
        variant_class: optional label (e.g. "indel") appended to the title so the
            SNV and indel views are distinguishable; empty value leaves the title as-is.
    """
    rcParams["font.size"] = 8
    x = sites.get_column("n_callers").fill_null(0).to_numpy()
    y = sites.get_column("cells_supporting_alt").to_numpy()

    # n_callers is a small integer; jitter x so co-located sites stay visible.
    jitter = (np.random.default_rng(0).random(len(x)) - 0.5) * 0.3

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(x + jitter, y, s=18, alpha=0.4, edgecolor="none", color="#3b6ea5")
    ax.axvline(ncallers_high - 0.5, color="grey", ls="--", lw=0.8)
    ax.axhline(multicell_min - 0.5, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("callers supporting site (n_callers)")
    ax.set_ylabel("cells supporting ALT")
    ax.set_xticks(range(0, int(x.max()) + 1))

    n_corner = int(((x >= ncallers_high) & (y >= multicell_min)).sum())
    class_tag = f" ({variant_class})" if variant_class else ""
    ax.set_title(
        f"Per-site per-cell vs abundance{class_tag}\n{sample_name} "
        f"(high-confidence corner: {n_corner}/{len(x)} sites)"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=120)
    plt.close(fig)


def _write_single_set_placeholder(
    output_path: str, suptitle: str, categories: set[str]
) -> None:
    """Write a labelled placeholder in place of an UpSet with fewer than two sets."""
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.axis("off")
    note = (
        f"single set: {', '.join(sorted(categories))}\nno intersection plot"
        if categories
        else "no sets present\nno intersection plot"
    )
    ax.text(0.5, 0.5, note, ha="center", va="center", fontsize=9)
    fig.suptitle(suptitle)
    fig.savefig(output_path, dpi=120)
    plt.close(fig)


def plot_generic_upset_plot(
    bool_count: pl.DataFrame,
    sample_name: str,
    title: str,
    modularity: str,
    discarded_members: List[str] = [],
    min_pct_show: int | None = 1,
    highlight_min_degree: int | None = None,
    highlight_label: str = "Usable",
    output_path: str | None = None,
    show_counts: bool = False,
    legend_loc: str = "best",
    legend_bbox_to_anchor: tuple[float, float] | None = None,
):
    """Plot an upset plot from a generic boolean counted dataframe.

    Args:
        bool_count (pl.DataFrame): Columns are flags, rows are true/false, the last column is the read count for the conditions in the row.
        sample_name (str): sample name
        title (str): Title of the plot except the sample name
        modularity (str): output-name prefix used when output_path is not given.
        discarded_members (List[str]): flags excluded from the highlight condition (still shown as categories).
        min_pct_show (int | None): minimum percentage of subset size to show; None disables the threshold.
        highlight_min_degree (int | None): when set, highlight intersections of at least this degree (e.g. NCALLERS>=2); otherwise highlight the "all-flags-present" intersection.
        highlight_label (str): label for the highlighted-subset legend entry.
        output_path (str | None): explicit output filename; defaults to "{modularity}_upset_plot.png".
        show_counts (bool): label bars with the absolute intersection size instead of the percentage. Defaults to False (percentages), keeping existing callers unchanged.
        legend_loc (str): matplotlib `loc` for the highlight legend; upsetplot's own "best" overlaps the bars.
        legend_bbox_to_anchor (tuple[float, float] | None): legend anchor in intersection-axes coordinates; negative x moves it out of the bar panel.
    """
    ids = bool_count.columns
    if READS in ids:
        ids.remove(READS)
    filtered_ids = [i for i in ids if i not in discarded_members]
    rcParams["font.size"] = 8
    memberships, read_counts = convert_df_to_upset_data(bool_count)
    memberships_hrn = convert_membership_names(memberships=memberships)
    hrn_ids = [label_for(i) for i in filtered_ids]
    final_output_path = output_path or f"{modularity}_upset_plot.png"

    # upsetplot can only build a MultiIndex from two or more distinct categories.
    # With <2 present (e.g. a sample with exactly one variant caller) UpSet raises
    # "'Index' object has no attribute 'levels'"; generate a placeholder instead
    # with a helper function
    distinct_categories = {c for group in memberships_hrn for c in group}
    if len(distinct_categories) < 2:
        _write_single_set_placeholder(
            final_output_path, f"{title}\n{sample_name}", distinct_categories
        )
        return

    upset_data = from_memberships(memberships=memberships_hrn, data=read_counts)
    # show_counts prints the absolute subset size on each bar; otherwise fall back
    # to the legacy percentage label so existing (barcode/amplicon) upsets are unchanged.
    upset_kwargs: dict = (
        {"show_counts": True} if show_counts else {"show_percentages": True}
    )
    if min_pct_show is not None:
        upset_kwargs["min_subset_size"] = f"{min_pct_show}%"
    upset = UpSet(upset_data, **upset_kwargs)
    if highlight_min_degree is not None:
        upset.style_subsets(
            min_degree=highlight_min_degree,
            facecolor="green",
            label=highlight_label,
        )
    else:
        upset.style_subsets(
            present=hrn_ids,
            facecolor="green",
            label=highlight_label,
        )
    axes = upset.plot()
    legend = axes["intersections"].get_legend()
    if legend is not None:
        legend.set_loc(legend_loc)
        if legend_bbox_to_anchor is not None:
            legend.set_bbox_to_anchor(legend_bbox_to_anchor)
    plt.suptitle(f"{title}\n{sample_name}")
    plt.savefig(final_output_path, dpi=120)
    plt.close()
