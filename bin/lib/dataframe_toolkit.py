import polars as pl
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix


def convert_polars_to_adata(
    df: pl.DataFrame,
    barcode_index_name: str,
    feature_index_name: str,
    values: str,
    annotations: pd.DataFrame | None = None,
) -> sc.AnnData:
    """Convert from a wide polars data frame to a scanpy adata object with index and vars.
    Also make the count matrix sparse csr_matrix

    """
    df = df.pivot(index=barcode_index_name, on=feature_index_name, values=values)
    df = df.fill_null(0)
    counts = csr_matrix(df.drop(barcode_index_name).to_numpy().astype("int32"))
    adata = sc.AnnData(counts)
    adata.obs_names = df[barcode_index_name].to_list()
    adata.var_names = df.columns[1:]
    if annotations is not None:
        adata.obs = adata.obs.join(annotations, how="left")
    return adata


def convert_adata_obs_to_polars(obs: pd.DataFrame):
    counts = obs.copy()[["total_counts"]]
    counts.reset_index(inplace=True, names=["barcode"])
    df = (
        pl.from_pandas(counts)
        .with_columns(pl.col("total_counts").cast(pl.UInt32))
        .rename({"total_counts": "rna_counts"})
    )
    return df


def get_most_freq_df(df: pl.LazyFrame, key) -> pl.LazyFrame:
    """Get the most frequent item of a column in a df

    Args:
        df (_type_): dataframe
        key (_type_): name of the column

    Returns:
        str: most frequent value
    """
    most_freq = df.group_by(key).agg(pl.len()).sort("len").tail(1).select(key)
    return most_freq
