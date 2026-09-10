import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_string_dtype,
)


def add_empty_cells(adata: ad.AnnData, new_cell_names: list[str]) -> ad.AnnData:
    """
    Add new cells with zero counts to an AnnData object.

    New obs metadata is filled with type-specific defaults:
    - integers -> 0
    - floats -> 0.0
    - booleans -> False
    - strings / objects / categoricals -> "na"

    Notes:
    - Existing categorical obs columns are extended with "na" if needed.
    - Duplicate cell names are rejected.
    - X is extended with sparse zero rows.
    - Every layer is extended with sparse zero rows of its own dtype so concat below
    keeps the size low.
    - uns is carried over from the original, a concat would otherwise drop it.
    """

    if not new_cell_names:
        return adata.copy()

    if len(new_cell_names) != len(set(new_cell_names)):
        raise ValueError("`new_cell_names` contains duplicates.")

    existing = set(adata.obs_names)
    overlapping = existing.intersection(new_cell_names)
    if overlapping:
        examples = sorted(overlapping)[:5]
        raise ValueError(
            f"`new_cell_names` overlaps with existing obs names. Examples: {examples}"
        )

    if adata.X is None:
        raise ValueError("Cannot add empty cells because `adata.X` is None.")

    adata = adata.copy()

    new_obs = pd.DataFrame(index=pd.Index(new_cell_names, name=adata.obs_names.name))

    for col in adata.obs.columns:
        series = adata.obs[col]
        dtype = series.dtype

        if isinstance(dtype, pd.CategoricalDtype):
            if "na" not in series.cat.categories:
                adata.obs[col] = series.cat.add_categories(["na"])

            new_obs[col] = pd.Series(
                pd.Categorical(
                    ["na"] * len(new_cell_names),
                    categories=adata.obs[col].cat.categories,
                    ordered=adata.obs[col].cat.ordered,
                ),
                index=new_cell_names,
            )

        elif is_integer_dtype(dtype):
            new_obs[col] = pd.Series(0, index=new_cell_names, dtype=dtype)

        elif is_float_dtype(dtype):
            new_obs[col] = pd.Series(0.0, index=new_cell_names, dtype=dtype)

        elif is_bool_dtype(dtype):
            new_obs[col] = pd.Series(False, index=new_cell_names, dtype=dtype)

        elif is_string_dtype(dtype):
            new_obs[col] = pd.Series("na", index=new_cell_names, dtype=dtype)

        else:
            # Safer fallback for object or unsupported extension dtypes.
            new_obs[col] = pd.Series("na", index=new_cell_names, dtype="object")

    zero_X = sp.csr_matrix(
        (len(new_cell_names), adata.n_vars),
        dtype=adata.X.dtype,
    )

    new_adata = ad.AnnData(
        X=zero_X,
        obs=new_obs,
        var=adata.var.copy(),
    )

    for key, layer in adata.layers.items():
        shape = (len(new_cell_names), adata.n_vars)
        new_adata.layers[key] = (
            sp.csr_matrix(shape, dtype=layer.dtype)
            if sp.issparse(layer)
            else np.zeros(shape, dtype=layer.dtype)
        )

    return ad.concat(
        [adata, new_adata],
        axis=0,
        join="outer",
        merge="same",
        uns_merge="first",
        fill_value=0,
    )


def add_metadata(
    mdata: mu.MuData, samplesheet: pd.DataFrame, sample_name: str
) -> mu.MuData:
    """
    Add metadata to a MuData object.

    Parameters
    ----------
    mdata : MuData
        The MuData object to which metadata will be added.
    samplesheet : pd.DataFrame
        A DataFrame containing metadata for the samples. The index should be the sample names.
    sample_name : str
        The name of the sample for which metadata will be added.

    Returns
    -------
    MuData
        A new MuData object with the added metadata.
    """

    if not samplesheet.index.is_unique:
        raise ValueError("samplesheet column 'sample' must be unique")

    if sample_name not in samplesheet.index:
        raise ValueError(f"Sample {sample_name!r} not found in samplesheet")

    row = samplesheet.loc[sample_name]

    mdata.obs["sample"] = pd.Categorical([sample_name] * mdata.n_obs)

    for col, value in row.items():
        if isinstance(value, str):
            mdata.obs[col] = pd.Categorical([value] * mdata.n_obs)
        else:
            mdata.obs[col] = value
    return mdata
