#!/usr/bin/env python3
"""Merge each per-cell pileup table with its genotype sibling.

    raw_variants.parquet       + variants_count.parquet       -> variants_snv.parquet
    raw_variants_indel.parquet + variants_indel_count.parquet -> variants_indel.parquet

Both sides describe the same cell-site observations, so this turns rows into
columns. Lossless: one row out per pileup row in, null genotype columns where
there was no call. reconciliation() is that statement as two counts.

Tiny minefields and ambigiuty to consider when working with this code:

- `filter` is read quality in the pileup table (lowBQ/lowMQ/lowDP) and call
  confidence in the genotype table (lowGQ). They disagree on most rows. Emitted
  as read_filter and gt_filter; never collapse them.
- a reference row is alt = ref in the pileup and alt IS NULL in the genotypes, so
  the REF join drops alt. Join on the full key and every hom-ref call silently
  vanishes. record_type carries the meaning instead.
- alt / alt_reads / vaf are null on REF rows: no ALT to describe. The pileup puts
  reference depth in them, which reads as a 100%-VAF carrier; ref_reads and
  total_reads say the same thing unambiguously.
- indel mean_bq exists only on the genotype side.
- genotypes are optional: count_variant*.py skips that file when it has nothing
  to call, and always writes the pileup one.
"""

import argparse
from typing import Final

import polars as pl
from lib.common_const import BARCODE, LibraryType, Modality
from lib.modality.dna.processing import (
    ALT_OUT,
    ALT_READS,
    CALLERS,
    CELLS_SUPPORTING_ALT,
    CELLS_TOTAL_AT_SITE,
    FEATURE,
    FILTER,
    GENOMIC_CHROM,
    GENOMIC_POS,
    GENOTYPE_LIKELIHOOD_SCHEMA,
    GQ,
    GT,
    IS_NO_CALL,
    MEAN_BQ,
    MEAN_MQ,
    N_CALLERS,
    OUTLIER_LOD,
    PL_HET,
    PL_HOM_ALT,
    PL_HOM_REF,
    POS_LOCAL,
    REF_OUT,
    REF_READS,
    TOTAL_READS,
    VAF,
    VARIANT_NAME,
    ZYGOSITY,
)
from lib.pipeline.status import ReasonCode, StatusRecord

OUT_SNV: Final[str] = "variants_snv.parquet"
OUT_INDEL: Final[str] = "variants_indel.parquet"

# Columns this merge introduces.
SAMPLE: Final[str] = "sample"
KIND: Final[str] = "kind"
RECORD_TYPE: Final[str] = "record_type"
READ_FILTER: Final[str] = "read_filter"
GT_FILTER: Final[str] = "gt_filter"
GENOTYPED: Final[str] = "genotyped"

KIND_SNV: Final[str] = "SNV"
KIND_INDEL: Final[str] = "INDEL"
RECORD_ALT: Final[str] = "ALT"
RECORD_REF: Final[str] = "REF"

# ALT rows join per allele; REF rows join without one.
KEY_ALT: Final[list[str]] = [BARCODE, FEATURE, POS_LOCAL, REF_OUT, ALT_OUT]
KEY_REF: Final[list[str]] = [BARCODE, FEATURE, POS_LOCAL]

# What the genotype table contributes, and the dtypes to fake when it is absent.
GENOTYPE_COLUMNS: Final[dict[str, pl.DataType]] = {
    VARIANT_NAME: pl.String,
    GT_FILTER: pl.String,
    **GENOTYPE_LIKELIHOOD_SCHEMA,
    ZYGOSITY: pl.String,
}
INDEL_GENOTYPE_COLUMNS: Final[dict[str, pl.DataType]] = {
    **GENOTYPE_COLUMNS,
    MEAN_BQ: pl.Float64,  # the indel pileup has none
}

# One schema for both outputs, so the two concatenate without a cast.
OUTPUT_COLUMNS: Final[list[str]] = [
    SAMPLE,
    BARCODE,
    KIND,
    RECORD_TYPE,
    FEATURE,
    POS_LOCAL,
    GENOMIC_CHROM,
    GENOMIC_POS,
    REF_OUT,
    ALT_OUT,
    VARIANT_NAME,
    TOTAL_READS,
    REF_READS,
    ALT_READS,
    VAF,
    MEAN_BQ,
    MEAN_MQ,
    READ_FILTER,
    GT_FILTER,
    IS_NO_CALL,
    CELLS_SUPPORTING_ALT,
    CELLS_TOTAL_AT_SITE,
    GT,
    GQ,
    PL_HOM_REF,
    PL_HET,
    PL_HOM_ALT,
    OUTLIER_LOD,
    ZYGOSITY,
    CALLERS,
    N_CALLERS,
    GENOTYPED,
]


def _attach(
    rows: pl.LazyFrame,
    genotypes: pl.LazyFrame | None,
    columns: dict[str, pl.DataType],
    on: list[str],
) -> pl.LazyFrame:
    """Left-join the genotype columns, or null them out if the sample has none."""
    if genotypes is None:
        return rows.with_columns(
            pl.lit(None, dtype=dtype).alias(name) for name, dtype in columns.items()
        )
    taken = [
        pl.col(FILTER).alias(GT_FILTER) if name == GT_FILTER else pl.col(name)
        for name in columns
    ]
    return rows.join(genotypes.select(*on, *taken), on=on, how="left")


def _shape(
    frame: pl.LazyFrame,
    sample: str,
    kind: str,
    record_type: str,
) -> pl.LazyFrame:
    """Label the rows and put them in the shared schema."""
    return frame.with_columns(
        pl.lit(sample, dtype=pl.String).alias(SAMPLE),
        pl.lit(kind, dtype=pl.String).alias(KIND),
        pl.lit(record_type, dtype=pl.String).alias(RECORD_TYPE),
        pl.col(GT).is_not_null().alias(GENOTYPED),
    ).select(OUTPUT_COLUMNS)


def merge_snv(
    raw: pl.LazyFrame,
    genotypes: pl.LazyFrame | None,
    sample: str,
) -> pl.LazyFrame:
    """raw_variants + variants_count. The ALT and REF halves join differently."""
    raw = raw.rename({FILTER: READ_FILTER})
    alt_rows = _attach(
        raw.filter(
            pl.col(ALT_OUT).is_not_null() & (pl.col(ALT_OUT) != pl.col(REF_OUT)),
        ),
        genotypes,
        GENOTYPE_COLUMNS,
        KEY_ALT,
    )
    ref_rows = _attach(
        raw.filter(pl.col(ALT_OUT) == pl.col(REF_OUT)),
        None if genotypes is None else genotypes.filter(pl.col(ALT_OUT).is_null()),
        GENOTYPE_COLUMNS,
        KEY_REF,
    ).with_columns(
        # callers, n_callers and the cell-support counts are already null here,
        # having joined on an allele no caller nominates; mean_bq / mean_mq stay,
        # being this row's own reads.
        pl.lit(None, dtype=pl.String).alias(ALT_OUT),
        pl.lit(None, dtype=pl.UInt32).alias(ALT_READS),
        pl.lit(None, dtype=pl.Float64).alias(VAF),
    )
    return pl.concat(
        [
            _shape(alt_rows, sample, KIND_SNV, RECORD_ALT),
            _shape(ref_rows, sample, KIND_SNV, RECORD_REF),
        ],
    )


def merge_indel(
    raw: pl.LazyFrame,
    genotypes: pl.LazyFrame | None,
    sample: str,
) -> pl.LazyFrame:
    """raw_variants_indel + variants_indel_count.

    No reference rows on either side: a cell that does not carry the event is a
    hom-ref genotype on the event's own row. So every row is ALT, and is_no_call
    is declared null to keep the two outputs on one schema.
    """
    raw = raw.rename({FILTER: READ_FILTER}).with_columns(
        pl.lit(None, dtype=pl.Boolean).alias(IS_NO_CALL),
    )
    rows = _attach(raw, genotypes, INDEL_GENOTYPE_COLUMNS, KEY_ALT)
    return _shape(rows, sample, KIND_INDEL, RECORD_ALT)


def reconciliation(
    merged: pl.DataFrame,
    raw_rows: int,
    genotype_rows: int,
) -> list[str]:
    """Losslessness as two counts. Empty means the merge held.

    Short of the pileup rows means a join key dropped some, over means one
    multiplied them. Short of the genotype rows means a call landed nowhere.
    """
    problems = []
    if merged.height != raw_rows:
        problems.append(
            f"merged rows {merged.height} != pileup rows in {raw_rows}: "
            "the join dropped or duplicated rows",
        )
    genotyped = int(merged[GENOTYPED].sum())
    if genotyped != genotype_rows:
        problems.append(
            f"genotyped rows {genotyped} != genotype rows in {genotype_rows}: "
            "a genotype matched no pileup row",
        )
    return problems


def _scan(path: str | None) -> pl.LazyFrame | None:
    return None if path is None else pl.scan_parquet(path)


def _rows(frame: pl.LazyFrame | None) -> int:
    return 0 if frame is None else int(frame.select(pl.len()).collect().item())


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge the per-cell pileup and genotype parquets.",
    )
    parser.add_argument("sample_name")
    parser.add_argument("snv_raw", help="raw_variants.parquet")
    parser.add_argument("indel_raw", help="raw_variants_indel.parquet")
    parser.add_argument("--snv-genotypes", help="variants_count.parquet, if any")
    parser.add_argument(
        "--indel-genotypes",
        help="variants_indel_count.parquet, if any",
    )
    args = parser.parse_args()

    status = StatusRecord(
        source_id=args.sample_name,
        library_type=LibraryType.DNA,
        modality=Modality.VARIANTS,
    )
    for out, merge, raw_path, genotype_path in (
        (OUT_SNV, merge_snv, args.snv_raw, args.snv_genotypes),
        (OUT_INDEL, merge_indel, args.indel_raw, args.indel_genotypes),
    ):
        raw, genotypes = _scan(raw_path), _scan(genotype_path)
        frame = merge(raw, genotypes, args.sample_name).collect()
        for problem in reconciliation(frame, _rows(raw), _rows(genotypes)):
            status.add(ReasonCode.INCONSISTENT_MERGE.value, f"{out}: {problem}")
        frame.write_parquet(out)
        print(f"{out}: {frame.height} rows, {int(frame[GENOTYPED].sum())} genotyped")

    if status.records_df.height:
        status.write_records("status.parquet")


if __name__ == "__main__":
    main()
