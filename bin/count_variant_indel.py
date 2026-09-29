#!/usr/bin/env python3
"""Per-cell attribution for INS/DEL calls — the indel sibling of count_variant.py.

Filtered to `type IN ('INS','DEL')`)

writes its own outputs next to the SNV ones:
  - raw_variants_indel.parquet  (per-cell per-allele support)
  - variants_indel_count.h5ad   (mod["indel"]: cells x called INS/DEL events)
  - variants_indel_count.parquet(per-cell genotype table)
"""

import argparse
from typing import Final

import polars as pl
import pysam
from lib.common_const import BARCODE, LibraryType, Modality
from lib.modality.dna.processing import (
    ALT_OUT,
    ALT_READS,
    CALLERS,
    CHROM,
    DEFAULT_ERROR_EPS,
    DEFAULT_ERROR_RHO,
    DEFAULT_FLAG_MIN_DP,
    DEFAULT_FLAG_MIN_MQ,
    DEFAULT_FLAG_MIN_VAF,
    DEFAULT_MIN_GQ,
    DEFAULT_MIN_SITE_POOLED_COVERAGE,
    FEATURE,
    FILTER,
    GENOMIC_CHROM,
    GENOMIC_POS,
    GQ,
    GT,
    MEAN_BQ,
    MEAN_MQ,
    N_CALLERS,
    OUTLIER_LOD,
    PL_HET,
    PL_HOM_ALT,
    PL_HOM_REF,
    POS,
    POS_LOCAL,
    REF_OUT,
    REF_READS,
    TOTAL_READS,
    VAF,
    VARIANT_NAME,
    ZYGOSITY,
    add_genomic_coords,
    build_read_cell_lut,
    empty_raw_indels,
    filter_indels,
    format_raw_indels_wide,
    indel_alt_long,
    indel_called_alleles,
    parse_amplicon_genomic_map,
    per_cell_indel_counts,
    per_site_indel_cell_support,
    unattributed_indel_site_rows,
)
from lib.modality.dna.variant_matrix import build_indel_anndata
from lib.pipeline.status import ReasonCode, StatusRecord

# Named per-cell thresholds are Nextflow params (conf/variant.config)
parser = argparse.ArgumentParser(
    description="Per-cell indel attribution + soft FILTER."
)
parser.add_argument("bam")
parser.add_argument("fasta")
parser.add_argument("sample_name")
parser.add_argument("barcode_mapping")
parser.add_argument("calls_parquet")
parser.add_argument("filtered_amplicon_reads")
parser.add_argument(
    "--min-site-pooled-coverage",
    type=int,
    default=DEFAULT_MIN_SITE_POOLED_COVERAGE,
    help="Minimum reads spanning an anchor, summed over all cells, to attribute it.",
)
parser.add_argument("--flag-min-dp", type=int, default=DEFAULT_FLAG_MIN_DP)
parser.add_argument("--flag-min-vaf", type=float, default=DEFAULT_FLAG_MIN_VAF)
# No --flag-min-bq here: an indel has no single anchor base quality to average from...
parser.add_argument("--flag-min-mq", type=float, default=DEFAULT_FLAG_MIN_MQ)
parser.add_argument(
    "--min-gq",
    type=float,
    default=DEFAULT_MIN_GQ,
    help="Genotype quality below which a per-cell call is labelled lowGQ / lowqual.",
)
parser.add_argument(
    "--error-eps",
    type=float,
    default=DEFAULT_ERROR_EPS,
    help="Per-base error rate used as the hom-ref / hom-alt expected ALT fraction.",
)
parser.add_argument(
    "--error-rho",
    type=float,
    default=DEFAULT_ERROR_RHO,
    help="Beta-binomial overdispersion (intra-class correlation); 0 is the binomial.",
)
args = parser.parse_args()

bam_path = args.bam
fasta_path = args.fasta
sample_name = args.sample_name
barcode_mapping_path = args.barcode_mapping
calls_parquet_path = args.calls_parquet
filtered_amplicon_reads_path = args.filtered_amplicon_reads

# Only INS and DEL variants are selected here
INDEL_TYPES: Final[list[str]] = ["INS", "DEL"]

# Every path that leaves without a matrix records why, so "this sample has no
# indel output" is answerable from status.csv rather than from the task log.
status = StatusRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.INDEL,
)


def _calls_lf() -> pl.LazyFrame:
    """Catalog parquet scoped to the indel records this script attributes."""
    return pl.scan_parquet(calls_parquet_path).filter(pl.col("type").is_in(INDEL_TYPES))


# Positions to pile up: the deduped union of every caller's indel position, lifted
# to the 0-based pileup coordinate. Identical selection to count_variant.py, only
# the type filter differs. No amplicon-anchor gate: on-target membership is
# enforced downstream by the filtered_amplicon_reads allowlist (start-OR-end
# anchored) + MIN_COVERAGE; gating on the start anchor only re-dropped 5'-edge /
# end-anchored-indel positions.
calls = _calls_lf().select(
    pl.col("chrom").alias(CHROM),
    pl.col("pos").cast(pl.Int64).alias(POS),
)
selected_features = (
    calls.with_columns(pl.col(POS) - 1)  # 0-based pileup column (samtools is 1-based)
    .rename({CHROM: FEATURE})
    .select(FEATURE, POS)
    .unique()  # dedupe positions shared by multiple callers; pileup once!
    .collect()
)

if selected_features.shape[0] == 0:
    # No indels called for this sample; still emit a schema-correct empty raw
    # output so the raw_variants_indel_parquet emit always has a file.
    empty_raw_indels().write_parquet("raw_variants_indel.parquet")
    status.record_and_exit(
        ReasonCode.LOW_FEATURES,
        f"No indels were called for {sample_name}, so there is nothing to "
        "attribute per cell.",
    )

samfile_handle = pysam.AlignmentFile(bam_path, "rb")
# Amplicon-to-genome lift from the self-describing FASTA headers (empty table for
# a plain FASTA, in which case genomic chrom/pos come out null).
# One major difference difference with count_variant for SNV is that no reffa handle
# needed: the indel identity is the read signature, not a reference base.
genomic_map = parse_amplicon_genomic_map(fasta_path)

# Map pileup reads to cells: read_name -> barcode, restricted to the reads that
# pass amplicon QC.
# This means peak memory follows the output (cells x signatures) rather than
# every (read, site) pair.
lut, cells = build_read_cell_lut(
    barcode_mapping_path=barcode_mapping_path,
    filtered_amplicon_reads_path=filtered_amplicon_reads_path,
)

print("Indel pileups started")
cell_counts = per_cell_indel_counts(
    selected_features=selected_features,
    samfile_handle=samfile_handle,
    lut=lut,
    cells=cells,
    min_site_pooled_coverage=args.min_site_pooled_coverage,
)
print("Indel pileups done")

if cell_counts is None or cell_counts.shape[0] == 0:
    empty_raw_indels().write_parquet("raw_variants_indel.parquet")
    status.record_and_exit(
        ReasonCode.LOW_COVERAGE,
        f"No called indel anchor in {sample_name} was spanned by a read from a "
        "called cell.",
    )

observed = cell_counts.lazy()
called_alleles = indel_called_alleles(_calls_lf())

# Per-cell support for each called allele, then per-site cell-support counts on
# the same (feature, pos_local, ref, alt) grouping the catalog joins on.
alt_long = indel_alt_long(observed=observed, called_alleles=called_alleles)
cell_support = per_site_indel_cell_support(
    observed=observed,
    alt_long=alt_long,
    flag_min_dp=args.flag_min_dp,
    flag_min_vaf=args.flag_min_vaf,
)

# Called indel alleles no cell carries a read of are surfaced here as
# zero-support site-level rows so nothing silently disappears and the drop
# stays countable.
attributed_wide = format_raw_indels_wide(
    alt_long=alt_long,
    cell_support=cell_support,
    flag_min_dp=args.flag_min_dp,
    flag_min_vaf=args.flag_min_vaf,
    flag_min_mq=args.flag_min_mq,
    genomic_map=genomic_map,
)
unattributed_wide = unattributed_indel_site_rows(
    called_alleles=called_alleles,
    alt_long=alt_long,
    genomic_map=genomic_map,
)
pl.concat([attributed_wide, unattributed_wide]).sink_parquet(
    "raw_variants_indel.parquet"
)

# No empty check: every piled-up anchor is a called one, so alt_long is empty
# only when cell_counts is, and that already exited above.
alt_long = alt_long.collect()

# Per-cell genotype per called allele from the beta-binomial likelihood. cells
# too shallow to separate the genotypes come out lowGQ / lowqual. Each indel
# allele is genotyped independently.
filtered_indels = filter_indels(
    alt_long=alt_long.lazy(),
    min_gq=args.min_gq,
    eps=args.error_eps,
    rho=args.error_rho,
).collect()

# mod["indel"]: cells x called INS/DEL events,
# X = AD/DP, the read evidence and the genotype in layers, and the event
# is added into structured var (event_type / start_pos / end_pos / ref_seq /
# alt_seq / event_len) instead of a VCF-anchored ref>alt string.
build_indel_anndata(
    observed=observed,
    alt_long=alt_long,
    genotypes=filtered_indels,
    called_alleles=called_alleles,
    calls=_calls_lf(),
).write_h5ad("variants_indel_count.h5ad")

# Per-cell genotype table with consensus support, lifted to genomic coords. Same
# schema as the SNV variants_count.parquet so the two union cleanly.
genotype = filtered_indels.lazy().join(
    called_alleles.select(FEATURE, POS_LOCAL, REF_OUT, ALT_OUT, CALLERS, N_CALLERS),
    on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT],
    how="left",
)
(
    add_genomic_coords(genotype, genomic_map)
    .select(
        BARCODE,
        FEATURE,
        POS_LOCAL,
        GENOMIC_CHROM,
        GENOMIC_POS,
        REF_OUT,
        ALT_OUT,
        TOTAL_READS,
        REF_READS,
        ALT_READS,
        VAF,
        MEAN_BQ,
        MEAN_MQ,
        GT,
        GQ,
        PL_HOM_REF,
        PL_HET,
        PL_HOM_ALT,
        OUTLIER_LOD,
        ZYGOSITY,
        VARIANT_NAME,
        FILTER,
        CALLERS,
        N_CALLERS,
    )
    .sink_parquet("variants_indel_count.parquet")
)
