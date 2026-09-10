#!/usr/bin/env python3
"""Caller-agreement QC for the per-sample variants parquet.

This script is a sibling of bin/qc_amplicon.py
It separates plotting from metric extraction so they can be independent of eachother.
Inputs match those of bin/variants_to_metrics_parquet.py, which is the long-format of parquet from VCF_TO_PARQUET

Produces:
  - {sample}_caller_upset.png: UpSet plot of variant caller per normalized SNV site;
    intersections of size >= 2 (NCALLERS >= 2) will be highlighted in green.
    Intersection bars are labelled with the absolute site count (not a percentage).
  - {sample}_caller_upset_indel.png: the same UpSet for INS/DEL sites (one combined
    indel plot). See the homopolymer caveat below.
  - metrics.parquet: flat MetricsRecord with `upset.<sorted_callers>`
    keys for SNVs and `upset_indel.<sorted_callers>` keys for indels (one per
    non-empty intersection, value = site count).
  - {sample}_persite_percell_vs_abundance.png: per-site scatter of caller agreement
    (n_callers) vs cell support (cells_supporting_alt) from the per-cell
    raw_variants.parquet (ADR-X2); the two-axis authenticity/abundance view.
  - {sample}_persite_percell_vs_abundance_indel.png: the same authenticity/abundance
    view for INS/DEL, driven by the per-cell raw_variants_indel.parquet.

None of these plots reach MultiQC

Indel caveat: the indel UpSet keys sites on (chrom, pos, ref, alt), so a single
biological indel represented at different lengths/anchors by different callers
(common in homopolymer/repeat runs) shows up as several caller-private sites.
The plot therefore over-states indel disagreement.

Usage:
    qc_variants.py <variants.parquet> <sample_name> <upset_png> <metrics.parquet> \\
        <raw_variants.parquet> <authenticity_png> \\
        <raw_variants_indel.parquet> <authenticity_indel_png> <indel_upset_png>
"""

import sys

import polars as pl
from lib.common_const import READS, LibraryType, Modality, label_for
from lib.metrics.metrics import MetricsRecord
from lib.plotting.common import plot_generic_upset_plot, plot_percell_variant_abundance

# TODO move this to a general config or fine tuning section
NCALLERS_HIGH_CONFIDENCE = 2

# Per-site key in the per-cell raw_variants.parquet.
SITE_KEY = ["feature", "pos_local_0based", "ref", "alt"]

# Catalog `type` values attributed to the (combined) indel UpSet.
INDEL_TYPES = ["INS", "DEL"]

# Legend placement for the caller UpSets, in intersection-axes coordinates.
UPSET_LEGEND_LOC = "upper right"
UPSET_LEGEND_BBOX_TO_ANCHOR = (-0.35, 1.0)


def render_authenticity_plot(
    raw_parquet_path: str,
    sample_name: str,
    output_path: str,
    variant_class: str = "",
) -> None:
    """Per-site per-cell vs abundance plot from a per-cell raw parquet.

    Works for both raw_variants.parquet (SNV) and raw_variants_indel.parquet
    (INS/DEL): both carry the same SITE_KEY plus n_callers / cells_supporting_alt
    columns, so only the title label (``variant_class``) differs.

    Deduplicates the cell-support columns to one row per site before
    plotting; writes an empty png when there is no per-cell support so
    the output channel stays shaped.
    """
    raw = pl.read_parquet(raw_parquet_path)
    if "cells_supporting_alt" not in raw.columns:
        open(output_path, "wb").close()
        return
    sites = (
        raw.filter(pl.col("cells_supporting_alt").is_not_null())
        .unique(subset=SITE_KEY)
        .select([*SITE_KEY, "n_callers", "cells_supporting_alt"])
    )
    if sites.is_empty():
        open(output_path, "wb").close()
        return
    plot_percell_variant_abundance(
        sites, sample_name, output_path, variant_class=variant_class
    )


def site_caller_booleans(df: pl.DataFrame, callers: list[str]) -> pl.DataFrame:
    """One row per (chrom, pos, ref, alt), one boolean column per caller."""
    sites = df.group_by(["chrom", "pos", "ref", "alt"]).agg(
        pl.col("caller").unique().alias("_callers")
    )
    sites = sites.with_columns(
        [pl.col("_callers").list.contains(c).alias(c) for c in callers]
    )
    return sites.select(callers)


def render_caller_upset(
    df: pl.DataFrame,
    sample_name: str,
    title: str,
    metric_prefix: str,
    record: MetricsRecord,
    output_path: str,
) -> None:
    """Caller-agreement UpSet for one (already type-filtered) catalog.

    Bars are labelled with the absolute intersection size (show_counts option),
    not a percentage. Generates one <metric_prefix>.<sorted_callers> metric per
    non-empty intersection. When the subset has no calls, writes an empty png
    and records nothing so the output channel stays shaped.
    """
    callers = sorted(df.get_column("caller").unique().to_list()) if df.height else []
    if not callers:
        open(output_path, "wb").close()
        return

    site_bools = site_caller_booleans(df, callers)
    bool_count = site_bools.group_by(callers).agg(pl.len().alias(READS))

    plot_generic_upset_plot(
        bool_count=bool_count,
        sample_name=sample_name,
        title=title,
        modularity="caller",
        min_pct_show=None,
        highlight_min_degree=NCALLERS_HIGH_CONFIDENCE,
        highlight_label=f"NCALLERS >= {NCALLERS_HIGH_CONFIDENCE}",
        output_path=output_path,
        show_counts=True,
        legend_loc=UPSET_LEGEND_LOC,
        legend_bbox_to_anchor=UPSET_LEGEND_BBOX_TO_ANCHOR,
    )

    for row in bool_count.iter_rows(named=True):
        present = sorted(c for c in callers if row[c])
        if not present:
            continue
        record.add(f"{metric_prefix}." + "_".join(present), row[READS])


def attribution_counts(
    catalog: pl.DataFrame, raw_path: str, types: list[str]
) -> tuple[int, int]:
    """(called, attributed) distinct-allele counts for one variant class.

    called = distinct alleles of types
    attributed = those that got a real per-cell row (non-null barcode) in the per-cell raw parquet.
    The variant catalog key is 1-based pos
    the per-cell key is 0-based pos_local_0based = pos - 1
    called - attributed is on TODO for silent-drop count: un-pileupablecalled
    alleles (indels also appear as zero-support null-barcode rows in the raw
    parquet. MNVs have no per-cell path, so attributed is 0 by construction).
    """
    called = (
        catalog.filter(pl.col("type").is_in(types))
        .select(
            pl.col("chrom").alias("feature"),
            (pl.col("pos") - 1).alias("pos_local_0based"),
            "ref",
            "alt",
        )
        .unique()
    )
    raw = pl.read_parquet(raw_path)
    if "barcode" in raw.columns:
        raw = raw.filter(pl.col("barcode").is_not_null())
    attributed = called.join(
        raw.select("feature", "pos_local_0based", "ref", "alt").unique(),
        on=["feature", "pos_local_0based", "ref", "alt"],
        how="inner",
    ).height
    return called.height, attributed


def main() -> None:
    parquet_path = sys.argv[1]
    sample_name = sys.argv[2]
    upset_png = sys.argv[3]
    out_yaml = sys.argv[4]
    raw_parquet_path = sys.argv[5]
    authenticity_png = sys.argv[6]
    indel_raw_parquet_path = sys.argv[7]
    authenticity_indel_png = sys.argv[8]
    indel_upset_png = sys.argv[9]

    catalog = pl.read_parquet(parquet_path)
    record = MetricsRecord(
        source_id=sample_name,
        library_type=LibraryType.DNA,
        modality=Modality.VARIANTS,
    )

    # per-site per-cell-vs-abundance plots are driven by the per-cell parquets and
    # are independent of caller presence, so render the SNV and indel views first.
    render_authenticity_plot(
        raw_parquet_path, sample_name, authenticity_png, variant_class="SNV"
    )
    render_authenticity_plot(
        indel_raw_parquet_path,
        sample_name,
        authenticity_indel_png,
        variant_class="indel",
    )

    # caller-agreement UpSets one per variant class.
    # Each handles its own no-calls case, so a sample with only SNVs (or only
    # indels) still gets an output for the other class.
    variants_label = label_for(Modality.VARIANTS.value)
    render_caller_upset(
        df=catalog.filter(pl.col("type") == "SNV"),
        sample_name=sample_name,
        title=f"{variants_label} caller agreement",
        metric_prefix="upset",
        record=record,
        output_path=upset_png,
    )
    render_caller_upset(
        df=catalog.filter(pl.col("type").is_in(INDEL_TYPES)),
        sample_name=sample_name,
        title=f"{variants_label} caller agreement (indel)",
        metric_prefix="upset_indel",
        record=record,
        output_path=indel_upset_png,
    )

    # Per-class attribution diagnostic: The goal is to findout how many called alleles reached
    # the per-cell layer. <class>_called - <class>_attributed is the silent-drop
    # count we'd like to track; MNV has no per-cell path so `mnvs_attributed` stays 0.
    for label, raw_path, types in (
        ("snvs", raw_parquet_path, ["SNV"]),
        ("indels", indel_raw_parquet_path, INDEL_TYPES),
        ("mnvs", indel_raw_parquet_path, ["MNV"]),
    ):
        called, attributed = attribution_counts(catalog, raw_path, types)
        record.add(f"{label}_called", called)
        record.add(f"{label}_attributed", attributed)

    record.write_records(out_yaml)


if __name__ == "__main__":
    main()
