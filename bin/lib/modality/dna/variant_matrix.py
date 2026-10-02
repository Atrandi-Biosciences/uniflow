"""Builds the two per-cell variant AnnData objects.

    mdata.mod["snv"]    cells x the amplicon positions that were piled up
    mdata.mod["indel"]  cells x called indel events

Both of these are projections of DFs count_variant.py and count_variant_indel.py
already produce. Nothing here changes or should change the pileup, the genotype
likelihood or the parquet outputs.

Every function is clearly defined what they do, the following is the general items
that is not immediately refelcted in each function:

  - the SNV axis is the measured positions, not every amplicon base: the pileup
    only visits called positions. uns["declared_target_sites"] keeps the rest.
  - DP is ACGT only. no-calls sit in AD_N, so DP is not raw_variants.total_reads
  - nonref_count is `DP - AD[ref_base]`, any nonref rather than the dominant one
  - GT is 1/2/3 for every genotyped pair whatever its GQ, -1 where there was
    nothing to genotype, and GQ is NaN exactly there. Filtering is the
    consumer's, with `GT in {1,2,3} AND GQ >= 30`.
  - GT and GQ are dense, every other layer sparse.
  - an indel del spans the deleted bases, so end - start + 1 is the deleted
    length.
  - an ins spans the two bases it sits between, so it is always 2
  - indel DP is anchor-locus depth, a superset of the reads spanning the
    deleted interval, so indel VAF is always conservative
"""

from collections.abc import Mapping
from typing import Final, NamedTuple

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pysam
import scipy.sparse as sp

from lib.common_const import BARCODE, LAYER_GQ, LAYER_GT, MODALITY, READS, Modality
from lib.modality.dna.processing import (
    ALT_OUT,
    ALT_READS,
    BASES,
    FEATURE,
    G0,
    GQ,
    GT,
    GT_HET,
    GT_HOM_ALT,
    GT_HOM_REF,
    INDEL_FEATURE,
    N_CALLERS,
    POS,
    POS_LOCAL,
    REF_BASE,
    REF_OUT,
    REF_READS,
    TARGET_GENOMIC_POS,
    TOTAL_READS,
    indel_anchor_coverage,
    parse_amplicon_headers,
)

########
# Global Variables and Mudata Grammar & Vocabulary
########

# The MuData slots two objects occupy.
SNV_MODALITY: Final[str] = Modality.SNV.value
INDEL_MODALITY: Final[str] = Modality.INDEL.value

# Layer names.
# AD_N is our invention. Standard ones are below
# This means no-calls stay visible outside DP.
LAYER_DP: Final[str] = "DP"
LAYER_AD: Final[str] = "AD"
LAYER_REF: Final[str] = "REF"
LAYER_OTHER_INDEL: Final[str] = "OTHER_INDEL"

# base to layer name map
# ACGT_LAYERS is the DP / nonref.
# N is outside of this
ACGT_LAYERS: Final[Mapping[str, str]] = {
    "A": "AD_A",
    "C": "AD_C",
    "G": "AD_G",
    "T": "AD_T",
}
NO_CALL_BASE: Final[str] = "N"
NO_CALL_LAYER: Final[str] = "AD_N"
SNV_AD_LAYERS: Final[Mapping[str, str]] = {**ACGT_LAYERS, NO_CALL_BASE: NO_CALL_LAYER}

# GT codes. -1 is "nothing to genotype", NOT "DP == 0": an indel cell whose
# every spanning read is anchor-deleted has DP > 0 and REF + AD == 0. 0 stays
# reserved so a stray sparse fill can never read as a genotype.
GT_NOT_GENOTYPED: Final[int] = -1
GT_HOM_REF_CODE: Final[int] = 1
GT_HET_CODE: Final[int] = 2
GT_HOM_ALT_CODE: Final[int] = 3
# Not 0: the model emits GQ 0.0000 on a zygosity boundary, so 0 would mean two
# things again.
GQ_NOT_GENOTYPED: Final[float] = np.nan

# var / obs field names.
VAR_AMPLICON_ID: Final[str] = "amplicon_id"
VAR_AMPLICON_POS: Final[str] = "amplicon_pos"
VAR_REF_BASE: Final[str] = "ref_base"
VAR_IS_MEASURED: Final[str] = "is_measured"
VAR_IS_TARGET_SITE: Final[str] = "is_target_site"
VAR_IS_CANDIDATE_VARIANT: Final[str] = "is_candidate_variant"
VAR_IS_CANDIDATE_INDEL: Final[str] = "is_candidate_indel"
VAR_N_CELLS_COVERED: Final[str] = "n_cells_covered"
VAR_N_CELLS_NONREF: Final[str] = "n_cells_with_nonref_support"
VAR_N_CELLS_INDEL: Final[str] = "n_cells_with_indel_support"
VAR_PB_DEPTH: Final[str] = "pseudobulk_depth"
VAR_PB_NONREF_COUNT: Final[str] = "pseudobulk_nonref_count"
VAR_PB_NONREF_VAF: Final[str] = "pseudobulk_nonref_vaf"
VAR_PB_MAJOR_NONREF: Final[str] = "pseudobulk_major_nonref"
VAR_EVENT_ID: Final[str] = "event_id"
VAR_EVENT_TYPE: Final[str] = "event_type"
VAR_START_POS: Final[str] = "start_pos"
VAR_END_POS: Final[str] = "end_pos"
VAR_REF_SEQ: Final[str] = "ref_seq"
VAR_ALT_SEQ: Final[str] = "alt_seq"
VAR_EVENT_LEN: Final[str] = "event_len"
VAR_EVENT_QC_PASS: Final[str] = "event_qc_pass"

OBS_TOTAL_DEPTH: Final[str] = "total_amplicon_depth"
OBS_N_SITES_COVERED: Final[str] = "n_sites_covered"
OBS_MEAN_SITE_DEPTH: Final[str] = "mean_site_depth"
OBS_MEDIAN_SITE_DEPTH: Final[str] = "median_site_depth"
OBS_N_SITES_NONREF: Final[str] = "n_sites_with_nonref_support"
OBS_FRACTION_SITES_NONREF: Final[str] = "fraction_sites_with_nonref_support"
OBS_N_EVENTS_COVERED: Final[str] = "n_indel_events_covered"
OBS_N_EVENTS_SUPPORTED: Final[str] = "n_indel_events_with_support"
OBS_TOTAL_INDEL_SUPPORT: Final[str] = "total_indel_support"
OBS_MEAN_EVENT_DEPTH: Final[str] = "mean_event_depth"
OBS_MEDIAN_EVENT_DEPTH: Final[str] = "median_event_depth"
OBS_FRACTION_EVENTS_SUPPORTED: Final[str] = "fraction_indel_events_with_support"

# Catalog of VCF_TO_PARQUET columns read for the pseudobulk var fields.
_CATALOG_CHROM: Final[str] = "chrom"
_CATALOG_POS: Final[str] = "pos"
_CATALOG_REF: Final[str] = "ref"
_CATALOG_ALT: Final[str] = "alt"
_CATALOG_CALLER: Final[str] = "caller"
_CATALOG_DP: Final[str] = "dp"
_CATALOG_AD_ALT: Final[str] = "ad_alt"

# Internal matrix-assembly columns.
_ROW: Final[str] = "_row"
_COL: Final[str] = "_col"
_NONREF: Final[str] = "_nonref"
_X: Final[str] = "_x"
_GT_CODE: Final[str] = "_gt_code"
# Reads at an anchor supporting anything other than the reference, i.e. every
# non-reference read before this column's own allele is taken out of it.
_ANCHOR_INDEL_READS: Final[str] = "_anchor_indel_reads"
_OTHER_INDEL: Final[str] = "_other_indel"
_VAR_NAME: Final[str] = "var_name"
_ABSENT_SEQ: Final[str] = "-"
# The allele identity an indel event joins on, everywhere.
_EVENT_KEYS: Final[list[str]] = [FEATURE, POS_LOCAL, REF_OUT, ALT_OUT]

# uns: the definitions that would otherwise have to be guessed.
# Added inside the object to make downstream analysis easier.
SNV_UNS: Final[dict[str, str]] = {
    "x_definition": "nonref_count / DP, nonref_count = DP - AD[ref_base]",
    "depth_definition": "AD_A + AD_C + AD_G + AD_T; no-calls are in AD_N",
    "measured_scope": (
        "the var axis is the positions at least one attributed read reached. "
        "The pileup only visits positions a caller called, so a position absent "
        "from var was never looked at, not looked at and found reference. Not "
        "the same as var['is_candidate_variant']: a called site can fail the "
        "pooled-coverage gate and end up unmeasured. "
        "uns['declared_target_sites'] carries the panel's declared targets and "
        "whether each was measured"
    ),
    "gt_encoding": (
        "-1 not genotyped, 1 hom_ref, 2 het, 3 hom_alt; 0 reserved, 4 reserved "
        "for multi-allelic and never emitted. Every genotyped pair gets 1/2/3 "
        "whatever its GQ, so a GQ-filtered view is "
        "`GT in (1,2,3) AND GQ >= min_gq`. -1 is not DP == 0: a cell whose "
        "every spanning read is anchor-deleted has depth but nothing to "
        "genotype. GQ is NaN exactly at -1, never 0, since 0 is a real GQ on a "
        "zygosity boundary"
    ),
    "obs_total_amplicon_depth": (
        "sum of per-position DP over the positions this modality measured, so "
        "read-observations rather than reads: a read spanning many measured "
        "positions counts once per position"
    ),
}
INDEL_UNS: Final[dict[str, str]] = {
    "x_definition": "AD / DP",
    "obs_total_amplicon_depth": (
        "sum of anchor DP over the distinct called anchors the cell covers, so "
        "read-observations rather than reads: a read spanning several called "
        "anchors counts once per anchor. Events sharing one anchor are counted "
        "once, not once per allele"
    ),
    "depth_definition": (
        "anchor-locus depth: every read spanning the anchor column, including "
        "reads whose anchor base is deleted by a larger event. A superset of "
        "the reads spanning the deleted interval, so VAF is conservative"
    ),
    "other_indel": (
        "reads at this anchor supporting neither the reference nor this "
        "column's event: a different indel signature (called or not), or a "
        "read whose anchor base is deleted by a larger event. AD + REF + "
        "OTHER_INDEL == DP exactly. It is limited to called anchors, since "
        "only those are piled up, so an un-called site has no column at "
        "all rather than a column of zeros"
    ),
    "gt_encoding": SNV_UNS["gt_encoding"],
    "event_axis": (
        "every called INS/DEL allele; a column of all zeros is a call no read "
        "could be attributed to (event_qc_pass False)"
    ),
    "cell_axis": (
        "every cell covering a called anchor, all of them genotyped: a covered "
        "non-carrier is hom_ref at alt_reads 0, not a blank. GT is -1 only "
        "where there was nothing to genotype, which includes the cell whose "
        "every spanning read is deleted through the anchor (DP > 0, "
        "REF + AD == 0)"
    ),
}

# Layer dtypes, kept as narrow as possible because GT and GQ are dense.
_COUNT_DTYPE: Final[str] = "int32"
_VAF_DTYPE: Final[str] = "float32"
_GT_DTYPE: Final[str] = "int8"
_GQ_DTYPE: Final[str] = "float32"


def gt_code_expr(gt_col: str = GT) -> pl.Expr:
    """Map the parquet genotype string to the layer code: 1/2/3, else -1.

    Blind to FILTER on purpose: GQ already carries the confidence. Read `gt`,
    never `zygosity`, which is nulled on a failed filter and would send every
    lowGQ row to -1.
    """
    return (
        pl.when(pl.col(gt_col) == GT_HOM_REF)
        .then(pl.lit(GT_HOM_REF_CODE))
        .when(pl.col(gt_col) == GT_HET)
        .then(pl.lit(GT_HET_CODE))
        .when(pl.col(gt_col) == GT_HOM_ALT)
        .then(pl.lit(GT_HOM_ALT_CODE))
        .otherwise(pl.lit(GT_NOT_GENOTYPED))
        .cast(pl.Int8)
        .alias(_GT_CODE)
    )


def _sparse_layer(
    frame: pl.DataFrame,
    value_col: str,
    shape: tuple[int, int],
    dtype: str,
) -> sp.csr_matrix:
    """Assemble one layer from the long (row, col, value) format.

    Zeros are dropped, so a stored entry always means a non-zero one. Read-count
    layers only; GT and GQ have no value free to serve as a fill.
    """
    matrix = sp.coo_matrix(
        (
            frame[value_col].to_numpy().astype(dtype),
            (frame[_ROW].to_numpy(), frame[_COL].to_numpy()),
        ),
        shape=shape,
        dtype=dtype,
    ).tocsr()
    matrix.eliminate_zeros()
    return matrix


def _dense_layer(
    frame: pl.DataFrame,
    value_col: str,
    shape: tuple[int, int],
    dtype: str,
    absent,
) -> np.ndarray:
    """Assemble GT / GQ as a dense array whose untouched slots read 'absent'."""
    matrix = np.full(shape, absent, dtype=dtype)
    matrix[frame[_ROW].to_numpy(), frame[_COL].to_numpy()] = (
        frame[value_col].to_numpy().astype(dtype)
    )
    return matrix


def _row_medians(matrix: sp.csr_matrix) -> np.ndarray:
    """Median of each row's stored (covered) values; 0 for an empty row."""
    medians = np.zeros(matrix.shape[0], dtype=_VAF_DTYPE)
    for row in range(matrix.shape[0]):
        start, stop = matrix.indptr[row], matrix.indptr[row + 1]
        if stop > start:
            medians[row] = np.median(matrix.data[start:stop])
    return medians


def _nnz(matrix: sp.csr_matrix, axis: int) -> np.ndarray:
    """Count covered entries: cells per column (axis 0), sites per row (axis 1)."""
    return matrix.getnnz(axis=axis).astype(_COUNT_DTYPE)


def _sum_per_row(matrix: sp.csr_matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel()


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Elementwise ratio, 0 where the denominator is 0 rather than NaN."""
    return np.divide(
        numerator,
        denominator,
        out=np.zeros(len(numerator), dtype=_VAF_DTYPE),
        where=denominator > 0,
    ).astype(_VAF_DTYPE)


class CellSummary(NamedTuple):
    """The per-cell statistics class for both modalities report, under neutral names.

    Each builder renames these to its own obs fields.
    """

    total_depth: np.ndarray
    n_covered: np.ndarray
    mean_depth: np.ndarray
    median_depth: np.ndarray
    n_signal: np.ndarray
    fraction_signal: np.ndarray


def _summarise_cells(depth: sp.csr_matrix, signal: sp.csr_matrix) -> CellSummary:
    covered = _nnz(depth, axis=1)
    total = _sum_per_row(depth)
    signalled = _nnz(signal, axis=1)
    return CellSummary(
        total_depth=total.astype(_COUNT_DTYPE),
        n_covered=covered,
        mean_depth=_safe_ratio(total, covered),
        median_depth=_row_medians(depth),
        n_signal=signalled,
        fraction_signal=_safe_ratio(signalled, covered),
    )


def _row_index(frame: pl.DataFrame) -> tuple[pd.Index, pl.DataFrame]:
    """Sorted barcode axis plus the barcode to row lookup map."""
    barcodes = frame[BARCODE].unique().sort()
    lookup = pl.DataFrame(
        {BARCODE: barcodes, _ROW: np.arange(len(barcodes), dtype="int64")}
    )
    return pd.Index(barcodes.to_list(), name=BARCODE), lookup


def _indexed_var(frame: pl.DataFrame) -> pd.DataFrame:
    """The var dataframe as pandas, indexed by var_name with no index label."""
    var = frame.to_pandas().set_index(_VAR_NAME)
    var.index.name = None
    return var


def _assemble(
    var: pd.DataFrame,
    x_matrix: sp.csr_matrix,
    layers: dict[str, sp.csr_matrix],
    obs: pd.DataFrame,
    modality: str,
    uns: dict,
) -> ad.AnnData:
    """Assemble the pieces into an AnnData with the pipeline's uns[MODALITY] key."""
    adata = ad.AnnData(X=x_matrix, obs=obs, var=var)
    for name, layer in layers.items():
        adata.layers[name] = layer
    adata.uns[MODALITY] = modality
    adata.uns[modality] = uns
    return adata


# ---------------------------------------------------------------------------
# SNV: cells x measured amplicon positions
# ---------------------------------------------------------------------------


def amplicon_position_axis(fasta_path: str) -> pl.DataFrame:
    """One row per position of every amplicon contig; _measured_axis cuts var from it.

    Carries amplicon_id, the 1-based amplicon_pos, ref_base, is_target_site and
    the column index. This is the only place the 0-based local position is
    converted to 1-based.
    """
    headers = parse_amplicon_headers(fasta_path)
    targets = (
        headers.filter(pl.col(TARGET_GENOMIC_POS).is_not_null())
        .select(
            pl.col(FEATURE).alias(VAR_AMPLICON_ID),
            (pl.col(TARGET_GENOMIC_POS) - pl.col(G0)).cast(pl.Int32).alias(POS_LOCAL),
            pl.lit(True).alias(VAR_IS_TARGET_SITE),
        )
        .unique()
    )

    with pysam.FastaFile(fasta_path) as handle:
        amplicon_ids: list[str] = []
        local_positions: list[int] = []
        ref_bases: list[str] = []
        for contig in handle.references:
            sequence = handle.fetch(contig).upper()
            amplicon_ids.extend([contig] * len(sequence))
            local_positions.extend(range(len(sequence)))
            ref_bases.extend(sequence)

    return (
        pl.DataFrame(
            {
                VAR_AMPLICON_ID: amplicon_ids,
                POS_LOCAL: local_positions,
                VAR_REF_BASE: ref_bases,
            },
            schema={
                VAR_AMPLICON_ID: pl.String,
                POS_LOCAL: pl.Int32,
                VAR_REF_BASE: pl.String,
            },
        )
        .join(targets, on=[VAR_AMPLICON_ID, POS_LOCAL], how="left")
        .with_columns(
            (pl.col(POS_LOCAL) + 1).cast(pl.Int32).alias(VAR_AMPLICON_POS),
            pl.col(VAR_IS_TARGET_SITE).fill_null(False),
        )
        .with_columns(
            pl.format("{}_{}", VAR_AMPLICON_ID, VAR_AMPLICON_POS).alias(_VAR_NAME),
            pl.int_range(pl.len(), dtype=pl.Int64).alias(_COL),
        )
    )


def snv_base_composition(raw_snp: pl.LazyFrame) -> pl.LazyFrame:
    """Per-cell base composition at the measured sites.

    One row per (barcode, feature, pos_local_0based) with AD_A/C/G/T/N, the
    ACGT-only depth, the non-reference count and the non-reference VAF.
    """
    # Both sides are uppercased here one from FASTA one from pileup.
    # On a soft-masked reference the two would disagree and every
    # cell would read as 100% non-reference.
    alleles = raw_snp.with_columns(
        pl.col(BASES).str.split(">").list.last().str.to_uppercase().alias(ALT_OUT),
        pl.col(REF_BASE).str.to_uppercase(),
    )
    composition = alleles.group_by(BARCODE, FEATURE, POS).agg(
        pl.col(REF_BASE).first(),
        *[
            pl.col(READS)
            .filter(pl.col(ALT_OUT) == base)
            .sum()
            .cast(pl.Int32)
            .alias(layer)
            for base, layer in ACGT_LAYERS.items()
        ],
        # Everything that is not a called ACGT base: N, and any other IUPAC code
        # the aligner produced. Catching them by exclusion rather than by name is
        # what keeps DP + AD_N equal to raw_variants.total_reads
        pl.col(READS)
        .filter(~pl.col(ALT_OUT).is_in(list(ACGT_LAYERS)))
        .sum()
        .cast(pl.Int32)
        .alias(NO_CALL_LAYER),
    )
    depth = pl.sum_horizontal([pl.col(layer) for layer in ACGT_LAYERS.values()])
    reference_depth = pl.sum_horizontal(
        [
            pl.when(pl.col(REF_BASE) == base).then(pl.col(layer)).otherwise(0)
            for base, layer in ACGT_LAYERS.items()
        ]
    )
    return (
        composition.with_columns(depth.cast(pl.Int32).alias(LAYER_DP))
        .with_columns(
            (pl.col(LAYER_DP) - reference_depth).cast(pl.Int32).alias(_NONREF)
        )
        .with_columns(
            pl.when(pl.col(LAYER_DP) > 0)
            .then(pl.col(_NONREF) / pl.col(LAYER_DP))
            .otherwise(0.0)
            .cast(pl.Float32)
            .alias(_X)
        )
        .rename({POS: POS_LOCAL})
    )


def snv_site_evidence(calls: pl.LazyFrame) -> pl.LazyFrame:
    """Per called position: n_callers plus the pseudobulk depth and support.

    Two consequences at a multi-allelic site, where the alleles can come from
    different callers over different denominators:
      - n_callers counts distinct callers at the POSITION, over all alleles.
        variants_count.parquet's identically named column is per allele, so the
        two can disagree and `n_callers >= 2` is not the same cut on each.
      - pseudobulk_nonref_vaf sums its numerator but takes a median denominator,
        so it is not bounded by 1. The per-cell X is, by construction.
    """
    # 1. Per allele. Several callers reported the same allele, each with its own
    # depth and alt-read count. Take the middle value of each
    per_allele = calls.group_by(
        pl.col(_CATALOG_CHROM).alias(VAR_AMPLICON_ID),
        (pl.col(_CATALOG_POS).cast(pl.Int32) - 1).alias(POS_LOCAL),
        pl.col(_CATALOG_ALT).alias(ALT_OUT),
    ).agg(
        pl.col(_CATALOG_DP).median().alias(VAR_PB_DEPTH),
        pl.col(_CATALOG_AD_ALT).median().alias(VAR_PB_NONREF_COUNT),
    )

    # 2. Per position. A position can carry more than one alternative allele.
    # Depth takes the middle of those; the alt counts get added, since they're
    # different alleles at the same spot
    per_site = per_allele.group_by(VAR_AMPLICON_ID, POS_LOCAL).agg(
        pl.col(VAR_PB_DEPTH).median(),
        pl.col(VAR_PB_NONREF_COUNT).sum(),
        pl.col(ALT_OUT)
        .sort_by(VAR_PB_NONREF_COUNT, descending=True, nulls_last=True)
        .first()
        .alias(VAR_PB_MAJOR_NONREF),
    )

    callers = calls.group_by(
        pl.col(_CATALOG_CHROM).alias(VAR_AMPLICON_ID),
        (pl.col(_CATALOG_POS).cast(pl.Int32) - 1).alias(POS_LOCAL),
    ).agg(pl.col(_CATALOG_CALLER).n_unique().cast(pl.Int32).alias(N_CALLERS))

    # Both sides group the same rows on the same keys, so every site has exactly
    # one row on each; a left join is enough and needs no coalescing.
    return (
        per_site.join(callers, on=[VAR_AMPLICON_ID, POS_LOCAL], how="left")
        .with_columns(
            pl.when(pl.col(VAR_PB_DEPTH) > 0)
            .then(pl.col(VAR_PB_NONREF_COUNT) / pl.col(VAR_PB_DEPTH))
            .cast(pl.Float32)
            .alias(VAR_PB_NONREF_VAF),
            pl.lit(True).alias(VAR_IS_CANDIDATE_VARIANT),
        )
        .select(
            VAR_AMPLICON_ID,
            POS_LOCAL,
            N_CALLERS,
            VAR_IS_CANDIDATE_VARIANT,
            VAR_PB_DEPTH,
            VAR_PB_NONREF_COUNT,
            VAR_PB_NONREF_VAF,
            VAR_PB_MAJOR_NONREF,
        )
    )


def _snv_entries(
    composition: pl.DataFrame,
    rows: pl.DataFrame,
    axis: pl.DataFrame,
    genotypes: pl.DataFrame,
    fasta_path: str,
) -> pl.DataFrame:
    """One row per (cell, measured position), carrying every layer's value.

    Raises if an attributed position is not in the FASTA, which would mean the
    calls and the reference disagree.
    """
    # Every join is validated via:
    # _sparse_layer builds a COO (coordinate) matrix, and COO sums duplicate
    # coordinates on tocsr(), so a single duplicated row here would
    # double a cell's counts and push GT into the unused code 4 rather
    # than raising an error.
    # The one-row-per-key guarantee lives entirely upstream, so it
    # is asserted rather than assumed (validate m:1 and 1:1).
    entries = (
        composition.join(rows, on=BARCODE, how="left", validate="m:1")
        .join(
            axis.select(
                pl.col(VAR_AMPLICON_ID).alias(FEATURE),
                pl.col(POS_LOCAL).cast(pl.Int64),
                _COL,
            ),
            on=[FEATURE, POS_LOCAL],
            how="left",
            validate="m:1",
        )
        .join(
            genotypes.select(
                BARCODE,
                FEATURE,
                pl.col(POS_LOCAL).cast(pl.Int64),
                gt_code_expr(),
                pl.col(GQ).cast(pl.Float32),
            ),
            on=[BARCODE, FEATURE, POS_LOCAL],
            how="left",
            validate="1:1",
        )
        .with_columns(
            pl.col(_GT_CODE).fill_null(GT_NOT_GENOTYPED),
            pl.col(GQ).fill_null(GQ_NOT_GENOTYPED),
        )
    )
    missing = entries.filter(pl.col(_COL).is_null())
    if missing.height:
        example = missing.select(FEATURE, POS_LOCAL).row(0)
        raise ValueError(
            f"{missing.height} attributed position(s) are absent from {fasta_path}; "
            f"first: feature={example[0]} pos_local_0based={example[1]}"
        )
    return entries


def _measured_axis(full_axis: pl.DataFrame, composition: pl.DataFrame) -> pl.DataFrame:
    """The var axis: the positions the pileup reached, re-indexed.

    The pileup visits called positions only (count_variant.py), so a dropped
    column was never looked at rather than looked at and found reference, and
    carried no measurement: is_candidate_variant False, n_callers 0, every
    pseudobulk_* null, zero cells covered. Its FASTA-derived fields survive in
    uns via _declared_target_sites.
    """
    return (
        full_axis.join(
            composition.select(
                pl.col(FEATURE).alias(VAR_AMPLICON_ID),
                pl.col(POS_LOCAL).cast(pl.Int32),
            ).unique(),
            on=[VAR_AMPLICON_ID, POS_LOCAL],
            how="semi",
        )
        .drop(_COL)
        .with_columns(pl.int_range(pl.len(), dtype=pl.Int64).alias(_COL))
    )


def _declared_target_sites(full_axis: pl.DataFrame, axis: pl.DataFrame) -> pd.DataFrame:
    """Every position the FASTA declares a target, and whether it was measured.

    var holds measured columns only, so without this an unmeasured target site
    leaves no trace in the object at all.
    """
    return (
        full_axis.filter(pl.col(VAR_IS_TARGET_SITE))
        .select(_VAR_NAME, VAR_AMPLICON_ID, VAR_AMPLICON_POS)
        .with_columns(pl.col(_VAR_NAME).is_in(axis[_VAR_NAME]).alias(VAR_IS_MEASURED))
        .to_pandas()
    )


def _snv_var(
    axis: pl.DataFrame,
    calls: pl.LazyFrame,
    depth: sp.csr_matrix,
    nonref: sp.csr_matrix,
) -> pd.DataFrame:
    """The SNV var: the measured axis, the catalog's site evidence, and the two
    per-column cell counts."""
    var = _indexed_var(
        axis.join(
            snv_site_evidence(calls).collect(),
            on=[VAR_AMPLICON_ID, POS_LOCAL],
            how="left",
        )
        .sort(_COL)
        .with_columns(
            pl.col(N_CALLERS).fill_null(0),
            pl.col(VAR_IS_CANDIDATE_VARIANT).fill_null(False),
            # Constant here, since the axis is the measured set. Kept for the
            # merged object, where measured in A not B is the whole question.
            pl.lit(True).alias(VAR_IS_MEASURED),
        )
        .select(
            _VAR_NAME,
            VAR_AMPLICON_ID,
            VAR_AMPLICON_POS,
            VAR_REF_BASE,
            VAR_IS_MEASURED,
            VAR_IS_TARGET_SITE,
            VAR_IS_CANDIDATE_VARIANT,
            N_CALLERS,
            VAR_PB_DEPTH,
            VAR_PB_NONREF_COUNT,
            VAR_PB_NONREF_VAF,
            VAR_PB_MAJOR_NONREF,
        )
    )
    var[VAR_N_CELLS_COVERED] = _nnz(depth, axis=0)
    var[VAR_N_CELLS_NONREF] = _nnz(nonref, axis=0)
    return var


def build_snv_anndata(
    raw_snp: pl.LazyFrame,
    genotypes: pl.DataFrame,
    calls: pl.LazyFrame,
    fasta_path: str,
) -> ad.AnnData:
    """The mod["snv"] object: cells x the measured amplicon positions.

    Args:
        raw_snp: the long per-allele pileup dataframe (format_raw_snps).
        genotypes: filter_snps output, one row per covered cell-site.
        calls: the catalog parquet filtered to type == 'SNV'.
        fasta_path: the amplicon reference; for the whole var axis.

    Columns are the positions the pileup reached; the declared targets it did
    not reach are in uns. Rows are the cells with at least one attributed read.
    MERGE_H5AD later widens them to the barcode union across modalities.
    """
    full_axis = amplicon_position_axis(fasta_path)
    composition = snv_base_composition(raw_snp).collect()
    axis = _measured_axis(full_axis, composition)
    obs_names, rows = _row_index(composition)
    entries = _snv_entries(composition, rows, axis, genotypes, fasta_path)

    shape = (len(obs_names), axis.height)
    layers = {
        **{
            layer: _sparse_layer(entries, layer, shape, _COUNT_DTYPE)
            for layer in SNV_AD_LAYERS.values()
        },
        LAYER_DP: _sparse_layer(entries, LAYER_DP, shape, _COUNT_DTYPE),
        LAYER_GT: _dense_layer(entries, _GT_CODE, shape, _GT_DTYPE, GT_NOT_GENOTYPED),
        LAYER_GQ: _dense_layer(entries, GQ, shape, _GQ_DTYPE, GQ_NOT_GENOTYPED),
    }
    nonref = _sparse_layer(entries, _NONREF, shape, _COUNT_DTYPE)
    depth = layers[LAYER_DP]

    cells = _summarise_cells(depth, signal=nonref)
    obs = pd.DataFrame(
        {
            OBS_TOTAL_DEPTH: cells.total_depth,
            OBS_N_SITES_COVERED: cells.n_covered,
            OBS_MEAN_SITE_DEPTH: cells.mean_depth,
            OBS_MEDIAN_SITE_DEPTH: cells.median_depth,
            OBS_N_SITES_NONREF: cells.n_signal,
            OBS_FRACTION_SITES_NONREF: cells.fraction_signal,
        },
        index=obs_names,
    )

    return _assemble(
        var=_snv_var(axis, calls, depth, nonref),
        x_matrix=_sparse_layer(entries, _X, shape, _VAF_DTYPE),
        layers=layers,
        obs=obs,
        modality=SNV_MODALITY,
        uns={
            **SNV_UNS,
            "declared_target_sites": _declared_target_sites(full_axis, axis),
        },
    )


# ---------------------------------------------------------------------------
# Indel: cells x called events
# ---------------------------------------------------------------------------


def _anchor_stripped(allele_col: str, other_len: pl.Expr) -> pl.Expr:
    """One allele with the shared VCF anchor removed; empty becomes "-"."""
    return pl.col(allele_col).str.slice(other_len).replace("", _ABSENT_SEQ)


def indel_event_var(called_alleles: pl.LazyFrame) -> pl.LazyFrame:
    """Turn VCF-anchored calls (AGG > A) into structured event fields.

    Adds event_type, start_pos / end_pos (1-based), ref_seq / alt_seq, the
    signed event_len, event_id and the unique var_name. Pure derivation from
    the catalog, INS / DEL only.

    A deletion occupies reference positions, so start..end are the deleted
    bases themselves and end - start + 1 is the deleted length. An insertion
    occupies none, so it is named by the two bases it sits between and
    end - start + 1 is always 2:

        AGG > A   del   ref_seq GG   alt_seq -    event_len -2
        A > AGG   ins   ref_seq -    alt_seq GG   event_len +2
    """
    len_ref = pl.col(REF_OUT).str.len_bytes().cast(pl.Int64)
    len_alt = pl.col(ALT_OUT).str.len_bytes().cast(pl.Int64)
    anchor = pl.col(POS_LOCAL).cast(pl.Int64) + 1
    is_insertion = len_alt > len_ref
    return (
        called_alleles.with_columns(
            pl.col(FEATURE).alias(VAR_AMPLICON_ID),
            pl.when(is_insertion)
            .then(pl.lit("ins"))
            .otherwise(pl.lit("del"))
            .alias(VAR_EVENT_TYPE),
            (len_alt - len_ref).cast(pl.Int32).alias(VAR_EVENT_LEN),
            _anchor_stripped(REF_OUT, len_alt).alias(VAR_REF_SEQ),
            _anchor_stripped(ALT_OUT, len_ref).alias(VAR_ALT_SEQ),
            pl.when(is_insertion)
            .then(anchor + len_ref - 1)
            .otherwise(anchor + len_alt)
            .cast(pl.Int32)
            .alias(VAR_START_POS),
            pl.when(is_insertion)
            .then(anchor + len_ref)
            .otherwise(anchor + len_ref - 1)
            .cast(pl.Int32)
            .alias(VAR_END_POS),
        )
        .with_columns(
            pl.format(
                "{}:{}:{}:{}:{}",
                VAR_AMPLICON_ID,
                VAR_EVENT_TYPE,
                VAR_START_POS,
                VAR_END_POS,
                VAR_ALT_SEQ,
            ).alias(VAR_EVENT_ID)
        )
        .with_columns(
            # var_names must be unique or AnnData slices wrongly. Two
            # anchors of one event share an event_id as well, and alleles
            # that are not pure indels (a complex REF>ALT, or GATK's `*`
            # spanning deletion) can additionally share the anchor, so the
            # anchor alone is not enough to separate them. The solution,
            # therefore, carries the whole allele key, which is unique by
            # construction: indel_called_alleles groups on exactly
            # (feature, pos_local_0based, ref, alt).
            pl.when(pl.col(VAR_EVENT_ID).is_duplicated())
            .then(pl.format("{}@{}:{}>{}", VAR_EVENT_ID, POS_LOCAL, REF_OUT, ALT_OUT))
            .otherwise(pl.col(VAR_EVENT_ID))
            .alias(_VAR_NAME)
        )
    )


def indel_site_evidence(calls: pl.LazyFrame) -> pl.LazyFrame:
    """Per called indel allele: the pseudobulk depth and support the callers saw.

    Same "median across callers" rule as the SNV side explained above, but
    per allele: two  events can share an anchor.
    """
    return (
        calls.group_by(
            pl.col(_CATALOG_CHROM).alias(FEATURE),
            (pl.col(_CATALOG_POS).cast(pl.Int64) - 1).alias(POS_LOCAL),
            pl.col(_CATALOG_REF).alias(REF_OUT),
            pl.col(_CATALOG_ALT).alias(ALT_OUT),
        )
        .agg(
            pl.col(_CATALOG_DP).median().alias(VAR_PB_DEPTH),
            pl.col(_CATALOG_AD_ALT).median().alias(VAR_PB_NONREF_COUNT),
        )
        .with_columns(
            pl.when(pl.col(VAR_PB_DEPTH) > 0)
            .then(pl.col(VAR_PB_NONREF_COUNT) / pl.col(VAR_PB_DEPTH))
            .cast(pl.Float32)
            .alias(VAR_PB_NONREF_VAF)
        )
    )


def _indel_entries(
    coverage: pl.DataFrame,
    rows: pl.DataFrame,
    events: pl.DataFrame,
    alt_long: pl.DataFrame,
    genotypes: pl.DataFrame,
) -> pl.DataFrame:
    """One row per (cell, event at an anchor that cell covers).

    Coverage drives the join, so non-carriers are present with their depth;
    support and genotype are joined on for the carriers among them.

    OTHER_INDEL is per column, not per anchor: it is the anchor's non-reference
    reads minus this event's own support, so at a shared anchor each event sees
    the others in it. Non-negative by construction, since alt_reads counts a
    subset of the reads that count as non-reference.
    """
    # The events join is deliberately many-to-many: several called alleles can
    # share one anchor, and each needs its own column for this cell. The others
    # must not fan out, and are validated because _sparse_layer's COO assembly
    # sums duplicate coordinates instead of raising on them.
    return (
        coverage.join(rows, on=BARCODE, how="left", validate="m:1")
        .join(events.select(*_EVENT_KEYS, _COL), on=[FEATURE, POS_LOCAL], how="inner")
        .join(
            alt_long.select(BARCODE, *_EVENT_KEYS, ALT_READS),
            on=[BARCODE, *_EVENT_KEYS],
            how="left",
            validate="m:1",
        )
        .join(
            # REF + AD == 0 at depth (every read another indel, or deleted
            # through the anchor) is filter_indels' flat-likelihood 0/0 at GQ 0,
            # not a genotype: left unjoined, it falls to -1 below.
            genotypes.filter(pl.col(REF_READS) + pl.col(ALT_READS) > 0).select(
                BARCODE,
                *_EVENT_KEYS,
                gt_code_expr(),
                pl.col(GQ).cast(pl.Float32),
            ),
            on=[BARCODE, *_EVENT_KEYS],
            how="left",
            validate="m:1",
        )
        .with_columns(
            pl.col(ALT_READS).fill_null(0),
            pl.col(_GT_CODE).fill_null(GT_NOT_GENOTYPED),
            pl.col(GQ).fill_null(GQ_NOT_GENOTYPED),
        )
        .with_columns(
            pl.when(pl.col(TOTAL_READS) > 0)
            .then(pl.col(ALT_READS) / pl.col(TOTAL_READS))
            .otherwise(0.0)
            .cast(pl.Float32)
            .alias(_X),
            (pl.col(_ANCHOR_INDEL_READS) - pl.col(ALT_READS)).alias(_OTHER_INDEL),
        )
    )


def build_indel_anndata(
    observed: pl.LazyFrame,
    alt_long: pl.DataFrame,
    genotypes: pl.DataFrame,
    called_alleles: pl.LazyFrame,
    calls: pl.LazyFrame,
) -> ad.AnnData:
    """The mod["indel"] object: cells x called indel events.

    Args:
        observed: per-cell indel signatures at the anchors (format_raw_indels).
        alt_long: the attributed (cell, called allele) support.
        genotypes: filter_indels output, carrying GT / GQ / FILTER.
        called_alleles: every called INS / DEL allele; this is the var axis.
        calls: the catalog parquet filtered to type IN ('INS','DEL').

    Columns are every called allele, so one nobody could attribute reads to is
    an all-zero column with event_qc_pass False. Rows are every cell covering a
    called anchor, all of them genotyped, so a covered non-carrier is a hom-ref
    with its depth visible. DP is anchor-locus depth, a superset of the reads
    spanning the deleted interval, so VAF is conservative.
    """
    events = (
        indel_event_var(called_alleles)
        .collect()
        .sort(_EVENT_KEYS)
        .with_columns(pl.int_range(pl.len(), dtype=pl.Int64).alias(_COL))
    )
    # AnnData only warns on duplicate var_names, then fails much later inside
    # ad.concat during the h5ad merge, with an error naming nothing indel
    # related. Fail here instead, where the offending allele is still in hand.
    collisions = events.filter(pl.col(_VAR_NAME).is_duplicated())
    if collisions.height:
        raise ValueError(
            f"{collisions.height} indel event(s) share a var_name: "
            f"{collisions[_VAR_NAME].unique().to_list()[:3]}"
        )

    coverage = (
        indel_anchor_coverage(observed)
        .with_columns(
            (pl.col(TOTAL_READS) - pl.col(REF_READS)).alias(_ANCHOR_INDEL_READS)
        )
        .collect()
    )
    obs_names, rows = _row_index(coverage)
    entries = _indel_entries(coverage, rows, events, alt_long, genotypes)

    shape = (len(obs_names), events.height)
    layers = {
        LAYER_AD: _sparse_layer(entries, ALT_READS, shape, _COUNT_DTYPE),
        LAYER_DP: _sparse_layer(entries, TOTAL_READS, shape, _COUNT_DTYPE),
        LAYER_REF: _sparse_layer(entries, REF_READS, shape, _COUNT_DTYPE),
        # Everything at the anchor that is neither reference nor this event, so
        # the three read layers partition DP. Still limited to called anchors:
        # an un-called cut site has no column here at all, rather than a zero.
        LAYER_OTHER_INDEL: _sparse_layer(entries, _OTHER_INDEL, shape, _COUNT_DTYPE),
        LAYER_GT: _dense_layer(entries, _GT_CODE, shape, _GT_DTYPE, GT_NOT_GENOTYPED),
        LAYER_GQ: _dense_layer(entries, GQ, shape, _GQ_DTYPE, GQ_NOT_GENOTYPED),
    }
    depth = layers[LAYER_DP]
    support = layers[LAYER_AD]

    return _assemble(
        var=_indel_var(events, calls, depth, support),
        x_matrix=_sparse_layer(entries, _X, shape, _VAF_DTYPE),
        layers=layers,
        obs=_indel_obs(
            depth,
            support,
            obs_names,
            _anchor_depth_per_cell(coverage, rows, len(obs_names)),
        ),
        modality=INDEL_MODALITY,
        uns=INDEL_UNS,
    )


def _indel_var(
    events: pl.DataFrame,
    calls: pl.LazyFrame,
    depth: sp.csr_matrix,
    support: sp.csr_matrix,
) -> pd.DataFrame:
    """The indel var: the normalised event, the catalog's evidence, and the
    per-column cell counts."""
    cells_with_support = _nnz(support, axis=0)
    var = _indexed_var(
        events.join(indel_site_evidence(calls).collect(), on=_EVENT_KEYS, how="left")
        .sort(_COL)
        .select(
            _VAR_NAME,
            VAR_EVENT_ID,
            VAR_AMPLICON_ID,
            VAR_EVENT_TYPE,
            VAR_START_POS,
            VAR_END_POS,
            VAR_REF_SEQ,
            VAR_ALT_SEQ,
            VAR_EVENT_LEN,
            pl.col(POS_LOCAL).cast(pl.Int32),
            REF_OUT,
            ALT_OUT,
            N_CALLERS,
            pl.lit(True).alias(VAR_IS_CANDIDATE_INDEL),
            VAR_PB_DEPTH,
            VAR_PB_NONREF_COUNT,
            VAR_PB_NONREF_VAF,
            pl.format("{}:{}:{}>{}", FEATURE, POS_LOCAL, REF_OUT, ALT_OUT).alias(
                INDEL_FEATURE
            ),
        )
    )
    cells_covered = _nnz(depth, axis=0)
    var[VAR_N_CELLS_COVERED] = cells_covered
    var[VAR_N_CELLS_INDEL] = cells_with_support
    # False = cells were sequenced here and none carried the signature, i.e. the
    # indels_called / indels_attributed gap. An anchor that was never piled up
    # at all (pooled coverage too low, or outside the amplicon allowlist) is a
    # different situation and stays null, so counting the False rows does not
    # over-report the gap.
    var[VAR_EVENT_QC_PASS] = pd.array(
        [
            (support > 0) if covered > 0 else pd.NA
            for covered, support in zip(cells_covered, cells_with_support)
        ],
        dtype="boolean",
    )
    # The indel anchor was piled up. Until now the only proxy was event_qc_pass null.
    var[VAR_IS_MEASURED] = cells_covered > 0
    # Declared but unknown: an indel target is an interval, the FASTA header
    # carries a point.
    var[VAR_IS_TARGET_SITE] = pd.array([pd.NA] * events.height, dtype="boolean")
    return var


def _anchor_depth_per_cell(
    coverage: pl.DataFrame, rows: pl.DataFrame, n_cells: int
) -> np.ndarray:
    """Reads per cell, summed over the distinct anchors it covers.

    Deliberately not the DP layer's row sum. Several called events can share one
    anchor, and the layer carries that anchor's depth once per event column, so
    summing across columns would count the same reads once per called allele.
    """
    totals = np.zeros(n_cells, dtype=_COUNT_DTYPE)
    per_cell = (
        coverage.join(rows, on=BARCODE, how="left")
        .group_by(_ROW)
        .agg(pl.col(TOTAL_READS).sum())
    )
    totals[per_cell[_ROW].to_numpy()] = per_cell[TOTAL_READS].to_numpy()
    return totals


def _indel_obs(
    depth: sp.csr_matrix,
    support: sp.csr_matrix,
    obs_names: pd.Index,
    total_depth: np.ndarray,
) -> pd.DataFrame:
    # mean / median stay per covered event, so they take the layer's row sums;
    # only the reported total has to be de-duplicated across shared anchors.
    cells = _summarise_cells(depth, signal=support)
    return pd.DataFrame(
        {
            OBS_TOTAL_DEPTH: total_depth,
            OBS_N_EVENTS_COVERED: cells.n_covered,
            OBS_N_EVENTS_SUPPORTED: cells.n_signal,
            OBS_TOTAL_INDEL_SUPPORT: _sum_per_row(support).astype(_COUNT_DTYPE),
            OBS_MEAN_EVENT_DEPTH: cells.mean_depth,
            OBS_MEDIAN_EVENT_DEPTH: cells.median_depth,
            OBS_FRACTION_EVENTS_SUPPORTED: cells.fraction_signal,
        },
        index=obs_names,
    )
