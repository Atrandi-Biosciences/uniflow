import numpy as np
import polars as pl
from lib.common_const import IS_CELL, N_FEATURES, TOTAL_COUNTS, CONFIDENCE_SCORE


def infer_otsu_threshold(values: np.ndarray, floor: int) -> int:
    """Infer threshold from 1D count values using numeric Otsu on log-counts."""
    vals = np.asarray(values, dtype=np.float64)
    if vals.size == 0:
        return floor
    if np.all(vals == vals[0]):
        return max(floor, int(vals[0]))

    log_vals = np.log10(vals + 1.0)
    counts, bin_edges = np.histogram(log_vals, bins=256)
    if np.count_nonzero(counts) <= 1:
        return max(floor, int(np.median(vals)))

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    mean1 = np.cumsum(counts * bin_centers) / np.clip(weight1, 1, None)
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / np.clip(weight2[::-1], 1, None))[
        ::-1
    ]
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    idx = int(np.argmax(variance12))
    threshold = np.rint(10.0 ** bin_centers[idx] - 1.0)
    return int(np.clip(threshold, floor, np.max(vals)))


def call_cells(
    barcode_metrics: pl.DataFrame,
    count_floor: int,
    feature_count_floor: int,
    min_confidence: float,
) -> pl.DataFrame:
    """Classify barcodes into real cell vs background using unsupervised thresholds."""
    count_threshold = infer_otsu_threshold(
        barcode_metrics.get_column(TOTAL_COUNTS).to_numpy(), floor=count_floor
    )
    feature_threshold = infer_otsu_threshold(
        barcode_metrics.get_column(N_FEATURES).to_numpy(),
        floor=feature_count_floor,
    )

    count_ratio = (pl.col(TOTAL_COUNTS) + 1.0) / (count_threshold + 1.0)
    feature_ratio = (pl.col(N_FEATURES) + 1.0) / (feature_threshold + 1.0)
    count_score = 1.0 / (1.0 + (-2.0 * count_ratio.log10()).exp())
    feature_score = 1.0 / (1.0 + (-2.0 * feature_ratio.log10()).exp())

    classes = barcode_metrics.with_columns(
        ((count_score + feature_score) / 2).alias(CONFIDENCE_SCORE),
    ).with_columns((pl.col(CONFIDENCE_SCORE) >= min_confidence).alias(IS_CELL))

    return classes
