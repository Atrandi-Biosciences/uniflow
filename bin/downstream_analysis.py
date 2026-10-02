#!/usr/bin/env python3
import mudata as mu
import sys
import numpy as np
import polars as pl

from lib.plotting.common import convert_df_to_upset_data

sample_name = sys.argv[1]
mudata_path = sys.argv[2]


# TODO: implement modalities as input to validate at runtime
mudata = mu.read_zarr(mudata_path)
mudata["gene_expression"].obs["total_counts"] = mudata["gene_expression"].X.sum(axis=1)
mudata["amplicon"].obs["total_counts"] = mudata["amplicon"].X.sum(axis=1)
mudata["gene_expression"].obs["n_genes_by_counts"] = (
    mudata["gene_expression"].X > 0
).sum(axis=1)
mudata["amplicon"].obs["n_genes_by_counts"] = (mudata["amplicon"].X > 0).sum(axis=1)

mudata.obs["is_cell"] = np.where(
    (mudata["amplicon"].obs["is_cell"]) & (mudata["gene_expression"].obs["is_cell"]),
    "double_positive",
    np.where(
        (mudata["amplicon"].obs["is_cell"])
        & (~mudata["gene_expression"].obs["is_cell"]),
        "amplicon_only",
        np.where(
            (~mudata["amplicon"].obs["is_cell"])
            & (mudata["gene_expression"].obs["is_cell"]),
            "gene_only",
            "double_negative",
        ),
    ),
)


cells_pl = (
    pl.from_pandas(mudata.obs.reset_index())
    .filter(pl.col("is_cell") != "double_negative")
    .select("index", "gene_expression:is_cell", "amplicon:is_cell")
    .with_columns(value=True)
    .group_by("gene_expression:is_cell", "amplicon:is_cell")
    .agg(pl.len())
)

membershipps = convert_df_to_upset_data(cells_pl)
print(membershipps)

# cells = from_memberships([["gene_expression_cell", "amplicon_cell"], ["gene_expression_cell"], ["amplicon_cell"]], data = [1683, 209, 37])
# upset = UpSet(cells, show_counts="%d", show_percentages=True, element_size=60)
# upset.style_subsets(present=["gene_expression_cell", "amplicon_cell"], facecolor="green", label="Kept")
# upset.plot()
# plt.rcParams["figure.figsize"] = (10,6)
# plt.suptitle("Cell calling by modality")
# plt.show()
