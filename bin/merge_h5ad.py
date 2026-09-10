#!/usr/bin/env python3
import sys
from pathlib import Path
import mudata as mu
import pandas as pd

from lib.common_const import (
    MODALITY,
)

from lib.anndata.utils import add_empty_cells

mu.set_options(pull_on_update=False)
sample_name: str = sys.argv[1]
sample_sheet_path: Path = Path(sys.argv[2])
experiment_id: str = sys.argv[3]
h5ad_path_list: list = sys.argv[4:]


h5ad_mapping = {}

all_barcodes = set()
h5ads = []
for h5ad_path in h5ad_path_list:
    adata = mu.read_h5ad(h5ad_path, mod=None)
    all_barcodes.update(adata.obs_names)
    h5ads.append(adata)

for adata in h5ads:
    modality = adata.uns[MODALITY]
    # TODO: This is currently creating "defaults" for obs columns without being clearly documented at metrics creation time. We need to fix this in the future.
    full_adata = add_empty_cells(adata, list(all_barcodes - set(adata.obs_names)))
    full_data_ordered = full_adata[sorted(all_barcodes), :]
    h5ad_mapping[modality] = full_data_ordered

merged = mu.MuData(h5ad_mapping)
samplesheet = pd.read_csv(sample_sheet_path)
merged.pull_obs()
merged.pull_var()

merged.obs["sample_name"] = pd.Categorical([sample_name] * merged.n_obs)
merged.obs["experiment_name"] = pd.Categorical([experiment_id] * merged.n_obs)
# TODO: fix add_metadata
# merged = add_metadata(merged, samplesheet, sample_name)

merged.write_h5mu("counts.h5mu")
merged.write_zarr(store="counts.zarr")
