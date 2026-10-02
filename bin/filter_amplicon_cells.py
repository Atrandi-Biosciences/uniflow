#!/usr/bin/env python3

import sys
import polars as pl
from typing import Final
import matplotlib.pyplot as plt
from lib.dataframe_toolkit import convert_polars_to_adata


from lib.common_const import (
    BARCODE,
    MODALITY,
    LibraryType,
    Modality,
    READS,
    CHROM,
    N_FEATURES,
    TOTAL_COUNTS,
    CONFIDENCE_SCORE,
    IS_CELL,
)

from lib.plotting.barcode import (
    experimental_knee_plot,
    plot_feature_count_scatter,
    plot_rank_knee,
    plot_histograms,
)
from lib.barcode.calling import call_cells
from lib.metrics.metrics import MetricsRecord
from lib.pipeline.status import ReasonCode, StatusRecord

filtered_counts_parquet: str = sys.argv[1]  # counts
full_read_parquet: str = sys.argv[2]  # Full read parquet path
sample_name: str = sys.argv[3]
force_cells = sys.argv[4]

if force_cells == "null":
    force_cells = None
else:
    force_cells = int(force_cells)


COUNT_THRESHOLD_DEFAULT_AMPLICON: Final[int] = 20
FEATURE_THRESHOLD_DEFAULT: Final[int] = 1
CONFIDENCE_THRESHOLD_DEFAULT: Final[float] = 0.7


metrics = MetricsRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.AMPLICON,
)

status = StatusRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.AMPLICON,
)

filtered_counts_lf = pl.scan_parquet(filtered_counts_parquet)
full_read_lf = pl.scan_parquet(full_read_parquet)


sorted_per_cell_total = (
    filtered_counts_lf.group_by(BARCODE).agg(pl.sum(READS)).sort(READS, descending=True)
).collect()


barcode_metrics = (
    filtered_counts_lf.group_by(BARCODE)
    .agg(
        pl.col(CHROM).n_unique().alias(N_FEATURES),
        pl.sum(READS).alias(TOTAL_COUNTS),
    )
    .collect()
)


barcode_scores = call_cells(
    barcode_metrics=barcode_metrics,
    count_floor=COUNT_THRESHOLD_DEFAULT_AMPLICON,
    feature_count_floor=FEATURE_THRESHOLD_DEFAULT,
    min_confidence=CONFIDENCE_THRESHOLD_DEFAULT,
)


if force_cells is not None:
    top_cells_stats = barcode_scores.sort(TOTAL_COUNTS, descending=True).head(
        force_cells
    )
else:
    top_cells_stats = barcode_scores.filter(pl.col(IS_CELL))

experimental_knee_plot(
    data=barcode_scores,
    sample_name=sample_name,
    ranking_col=TOTAL_COUNTS,
    modality=Modality.AMPLICON.value,
)

plot_feature_count_scatter(
    data=barcode_scores,
    sample_name=sample_name,
    confidence_threshold=CONFIDENCE_THRESHOLD_DEFAULT,
)

top_cells = top_cells_stats.select(BARCODE)
metrics.add("filtered_cells", top_cells.shape[0])
total_reads = full_read_lf.select(pl.len()).collect().item()


filtered_reads_in_cells = (
    filtered_counts_lf.join(top_cells.lazy(), on=BARCODE)
    .select(pl.sum(READS))
    .collect()
    .item()
)
metrics.add("total_usable_reads_in_cells", filtered_reads_in_cells)

fraction_usable_reads_in_cells = (
    filtered_reads_in_cells / total_reads if total_reads > 0 else 0
)
metrics.add("fraction_usable_reads_in_cells", fraction_usable_reads_in_cells)


# TODO: Split the rest of this code into another process
if top_cells.shape[0] == 0:
    status.add(
        ReasonCode.NO_CELLS.value,
        "No cells were detected in the sample. Check the anchors and filtering criteria.",
    )
    status.write_records("status.parquet")
else:
    top_cells.write_parquet(f"{sample_name}_top_cells.parquet")
top_barcodes = top_cells.to_series().to_list()
plt.figure(figsize=(6, 4), dpi=120)
umi_threshold = plot_rank_knee(
    sorted_per_cell_total,
    knee_ncells=top_cells.shape[0],
    sample_name=sample_name,
    value_type=READS,
)

plot_histograms(
    sorted_per_cell_total,
    sample_name=sample_name,
    umi_threshold=umi_threshold,
    data_type=READS,
)

metrics.write_records("metrics.parquet")
n_rows = filtered_counts_lf.select(pl.len()).collect().item()

if n_rows < 20:
    print(
        "Warning: Less than 20 barcodes with reads after filtering. Check the anchors and filtering criteria."
    )
    sys.exit(0)

adata = convert_polars_to_adata(
    df=filtered_counts_lf.collect(),
    barcode_index_name=BARCODE,
    feature_index_name=CHROM,
    annotations=barcode_scores.select(BARCODE, CONFIDENCE_SCORE, IS_CELL)
    .to_pandas()
    .set_index(BARCODE),
    values=READS,
)
adata.uns[MODALITY] = Modality.AMPLICON.value

# adata.obs[IS_CELL] = pd.Series(
#     [idx in top_barcodes for idx in adata.obs.index],
#     index=adata.obs.index,
#     name=IS_CELL,
# ).astype(bool)

adata.write_h5ad(f"filtered_{Modality.AMPLICON.value}.h5ad")
