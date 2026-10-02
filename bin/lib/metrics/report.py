import re
from pathlib import Path

import polars as pl


def slug(s: str) -> str:
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


def build_config(
    metrics_reference_csv: Path,
    aggregated_metrics_csv: Path,
    scope: str,
    internal_flag: bool,
) -> dict:
    metrics = pl.read_csv(aggregated_metrics_csv)
    # Back-compat: old modality name
    metrics = metrics.with_columns(
        pl.when(pl.col("modality") == "none")
        .then(pl.lit("barcode"))
        .otherwise(pl.col("modality"))
        .alias("modality")
    )

    ref = pl.read_csv(
        metrics_reference_csv,
        schema_overrides={
            "library_type": pl.Utf8,
            "modality": pl.Utf8,
            "metric_key": pl.Utf8,
            "human_readable_name": pl.Utf8,
            "description": pl.Utf8,
            "report": pl.Utf8,
            "customer_facing": pl.Boolean,
            "min_warning_threshold": pl.Utf8,
            "max_warning_threshold": pl.Utf8,
            "min_failure_threshold": pl.Utf8,
            "max_failure_threshold": pl.Utf8,
        },
    )
    # if customer facing, delete internal metrics
    if not internal_flag:
        ref = ref.filter(pl.col("customer_facing"))

    ref = ref.with_columns(pl.col("report").str.to_lowercase())

    if scope not in {"aggregated", "single"}:
        raise ValueError("scope must be 'aggregated' or 'single'")

    # Include rows where report == scope or report == both
    ref = ref.filter(pl.col("report").is_in([scope, "both"]))

    # Stable order from source file
    ref = ref.with_row_index("ref_order")

    merged = metrics.join(
        ref,
        on=["library_type", "modality", "metric_key"],
        how="inner",
    )

    sections = {}
    order = []

    # Unique combos ordered by first appearance in reference
    ref_combos = (
        ref.select(["library_type", "modality", "ref_order"])
        .group_by(["library_type", "modality"])
        .agg(pl.col("ref_order").min().alias("ref_order"))
        .sort("ref_order")
    )

    for combo in ref_combos.iter_rows(named=True):
        lib = combo["library_type"]
        mod = combo["modality"]

        sub = merged.filter(
            (pl.col("library_type") == lib) & (pl.col("modality") == mod)
        )
        if sub.height == 0:
            continue

        ref_cols = ref.filter(
            (pl.col("library_type") == lib) & (pl.col("modality") == mod)
        ).sort("ref_order")

        section_id = f"metrics_{slug(lib)}_{slug(mod)}"
        metric_keys = ref_cols.get_column("metric_key").to_list()

        wide = sub.pivot(
            values="value",
            index="source_id",
            on="metric_key",
            aggregate_function="first",
        ).rename({"source_id": "Sample"})

        cols_present = ["Sample"] + [k for k in metric_keys if k in wide.columns]
        wide = wide.select(cols_present)

        wide.write_csv(f"{section_id}.tsv", separator="\t")
        print(f"wrote out {section_id}")

        headers = {"Sample": {"title": "Sample"}}

        for r in ref_cols.iter_rows(named=True):
            key = r["metric_key"]
            if key not in wide.columns:
                continue

            header = {
                "title": r["human_readable_name"],
                "description": r["description"],
                "format": (
                    r["format"]
                    if r["format"]
                    else ("{:.3f}" if key.startswith("fraction_") else "{:,.0f}")
                ),
            }

            warn_rules = []
            fail_rules = []
            if r["min_warning_threshold"]:
                warn_rules.append({"lt": float(r["min_warning_threshold"])})
            if r["min_failure_threshold"]:
                fail_rules.append({"lt": float(r["min_failure_threshold"])})
            if r["max_warning_threshold"]:
                warn_rules.append({"gt": float(r["max_warning_threshold"])})
            if r["max_failure_threshold"]:
                fail_rules.append({"gt": float(r["max_failure_threshold"])})

            if warn_rules or fail_rules:
                header["cond_formatting_rules"] = {}
                if warn_rules:
                    header["cond_formatting_rules"]["warn"] = warn_rules
                if fail_rules:
                    header["cond_formatting_rules"]["fail"] = fail_rules

            headers[key] = header

        mod_title = mod.replace("_", " ").title()
        title = f"{lib} {mod_title} Metrics"

        sections[section_id] = {
            "file_format": "tsv",
            "section_name": title,
            "description": f"Aggregated {lib} {mod.replace('_', ' ')} metrics per sample or library.",
            "plot_type": "table",
            "parent_id": f"metrics_{slug(lib)}",
            "parent_name": f"{lib} Metrics",
            "parent_description": f"Custom metrics for {lib} libraries.",
            "pconfig": {
                "id": section_id,
                "title": title,
            },
            "headers": headers,
        }

        order.append(section_id)

    config = {
        "custom_content": {"order": order},
        "custom_data": sections,
        "sp": {sec_id: {"fn": f"{sec_id}.tsv"} for sec_id in sections.keys()},
    }

    return config
