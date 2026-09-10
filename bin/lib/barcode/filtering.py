import numpy as np
import polars as pl
from typing import Tuple
from lib.common_const import COUNTS, BARCODE


def knee_point(y: np.ndarray, smooth_window: int = 101) -> int:
    """
    Detect the "knee" in a monotonically decreasing curve.

    The knee is taken as the point of maximum curvature (sharp bend) in the
    smoothed trend.

    Parameters
    ----------
    y : np.ndarray
        1-D array of y-coordinates.  It is assumed the x-coordinates are
        simply the indices 0, 1, 2, … len(y)-1.
    smooth_window : int, optional (default=101)
        Size of the moving-average window used to smooth the curve before
        curvature computation.  The window is automatically truncated to the
        largest odd integer ≤ len(y) and never allowed to drop below 5.

    Returns
    -------
    int
        Index of the detected knee.  The value is clipped so that it lies
        between the 1st and 90th percentile positions of the curve, and
        always within the array bounds.

    Notes
    -----
    - The curve is padded at the edges by mirroring neighbouring points to
      avoid boundary artefacts during convolution.
    - Curvature is computed as |y''| / (1 + y'²)^(3/2), where derivatives
      are estimated by `np.gradient`.
    - If the input has fewer than 5 points, the last valid index is returned
      immediately.
    """
    n = len(y)

    # Too little data → return last valid index
    if n < 5:
        return n - 1

    # Ensure the smoothing window is odd and within reasonable bounds
    w = min(smooth_window, n if n % 2 == 1 else n - 1)
    if w < 5:
        w = 5 if n >= 5 else n
    half = w // 2

    # Mirror the head and tail of the array to avoid edge artefacts
    pad_left = y[1 : half + 1][::-1] if half > 0 else np.array([], dtype=y.dtype)
    pad_right = y[-half - 1 : -1][::-1] if half > 0 else np.array([], dtype=y.dtype)
    y_pad = np.concatenate([pad_left, y, pad_right])

    # Smooth the curve with a simple moving average
    y_smooth = np.convolve(y_pad, np.ones(w) / w, mode="valid")

    # First and second derivatives
    dy = np.gradient(y_smooth)
    d2y = np.gradient(dy)

    # Curvature formula: |y''| / (1 + y'²)^(3/2)
    curvature = np.abs(d2y) / np.power(1.0 + dy**2, 1.5)

    # Restrict the search to the central 1%–90% region to ignore noisy ends
    lo = int(0.01 * n)
    hi = max(lo + 1, int(0.90 * n))

    # Return the index of the maximum curvature, clipped to valid range
    knee_idx = lo + int(np.argmax(curvature[lo:hi]))
    return int(np.clip(knee_idx, 0, n - 1))


def knee_cells_from_umis(
    df: pl.DataFrame, expected_cells: int, min_umi: int = 100
) -> Tuple[pl.DataFrame, bool]:
    """Find the top cells based on a knee heuristic

    Args:
        df (pl.DataFrame): Dataframe of counts and barcode
        min_umi (int, optional): Min number of umi counts before finding the knee. Defaults to 100.

    Returns:
        pl.DataFrame | None: _description_
    """

    y = df.select(COUNTS).to_numpy()
    # y = np.sort(total_counts.astype(float))[::-1]
    y = y[y >= min_umi]
    if y.size == 0 or y.size <= (expected_cells / 5):
        print("Number of cells is too low, will use expected cells instead.")
        return df.select(BARCODE).head(expected_cells), False
    k = knee_point(y)

    return df.select(BARCODE).head(k + 1), True
