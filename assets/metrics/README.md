# Metrics Reference

## Purpose
`metrics_reference.csv` is the ground-truth dictionary for metric metadata used in the aggregated metrics report. It maps each metric key to a human-readable display name, a short description, and reporting/threshold behavior. This file is the authoritative source for metric labeling, descriptions, and thresholding. Update it whenever a metric is added, renamed, or removed.

## File Contents
Each row defines a single metric for a given `library_type` and `modality`.

Columns:
- `library_type`: The library type the metric applies to (e.g., `RNA`, `DNA`).
- `modality`: The modality context for the metric (e.g., `barcode`, `amplicon`, `gene_expression`).
- `metric_key`: The canonical metric identifier used in pipeline outputs.
- `human_readable_name`: Display name to use in reports.
- `description`: Concise explanation of the metric.
- `report`: A **flag** indicating where the metric should appear:
  - `aggregated` = show only in aggregated (multi-sample/library) report.
  - `single` = show only in single-sample/library report.
  - `both` = show in both aggregated and single report.
  - empty = do not include in any report.
- `min_warning_threshold`: Values **below** this trigger a warning.
- `max_warning_threshold`: Values **above** this trigger a warning.
- `min_failure_threshold`: Values **below** this trigger a failure.
- `max_failure_threshold`: Values **above** this trigger a failure.

