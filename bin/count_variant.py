#!/usr/bin/env python3
import argparse

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
    DEFAULT_FLAG_MIN_BQ,
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
    empty_raw_variants,
    filter_snps,
    format_raw_snps_wide,
    parse_amplicon_genomic_map,
    per_cell_snv_counts,
    per_site_cell_support,
)
from lib.modality.dna.variant_matrix import build_snv_anndata
from lib.pipeline.status import ReasonCode, StatusRecord

# Named per-cell thresholds are Nextflow params (conf/variant.config); the
# defaults here only apply when the script is run by hand.
parser = argparse.ArgumentParser(description="Per-cell SNV attribution + soft FILTER.")
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
    help="Minimum reads spanning a site, summed over all cells, to attribute it.",
)
parser.add_argument("--flag-min-dp", type=int, default=DEFAULT_FLAG_MIN_DP)
parser.add_argument("--flag-min-vaf", type=float, default=DEFAULT_FLAG_MIN_VAF)
parser.add_argument("--flag-min-bq", type=float, default=DEFAULT_FLAG_MIN_BQ)
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

status = StatusRecord(
    source_id=sample_name,
    library_type=LibraryType.DNA,
    modality=Modality.SNV,
)


# the per-cell pileup is caller-independent. limit the scope to SNV-only
# so indel/MNV rows in the catalog parquet are excluded here.
calls = (
    pl.scan_parquet(calls_parquet_path)
    .filter(pl.col("type") == "SNV")
    .select(
        pl.col("chrom").alias(CHROM),
        pl.col("pos").cast(pl.Int32).alias(POS),
    )
)

# Positions to pile up = the deduped union of every caller's called SNV positions,
# lifted to the 0-based pileup coordinate. No amplicon-anchor gate here: on-target
# reads is already enforced downstream by the filtered_amplicon_reads
# allowlist (itself start-OR-end anchored) and by MIN_COVERAGE, so gating on the
# start anchor only re-dropped 5'-edge / end-anchored-indel positions.
selected_features = (
    calls.with_columns(
        pl.col(POS)
        - 1,  # our pileup uses 0-based positions whereas samtools/VCF are 1-based
    )
    .rename({CHROM: FEATURE})
    .select(FEATURE, POS)
    .unique()  # dedupe positions shared by multiple callers; pileup once!
    .collect()
)
if selected_features.shape[0] > 0:
    samfile_handle = pysam.AlignmentFile(bam_path, "rb")
    ref_handle = pysam.FastaFile(fasta_path)
    # amplicon to genome lift happens using the self-describing FASTA headers
    # (pysam only reports contig IDs only, not the metadata description). Empty for a
    # plain FASTA, in which case the genomic chrom/pos columns will come out as "null".
    genomic_map = parse_amplicon_genomic_map(fasta_path)

    # Map pileup reads to cells: read_name -> barcode, restricted to the reads that
    # pass amplicon QC.
    # This means peak memory follows the output (cells x signatures) rather than
    # every (read, site) pair.
    lut, cells = build_read_cell_lut(
        barcode_mapping_path=barcode_mapping_path,
        filtered_amplicon_reads_path=filtered_amplicon_reads_path,
    )

    # For now, we filter based on selected bulk features. This is where the filtering happens
    print("Pileups started")
    cell_counts = per_cell_snv_counts(
        selected_features=selected_features,
        samfile_handle=samfile_handle,
        ref_handle=ref_handle,
        lut=lut,
        cells=cells,
        min_site_pooled_coverage=args.min_site_pooled_coverage,
    )
    print("Pileups done")

    if cell_counts is None or cell_counts.shape[0] == 0:
        # Every selected site failed the pooled coverage or none of the
        # spanning reads belongs to a called cell
        empty_raw_variants().write_parquet("raw_variants.parquet")
        status.record_and_exit(
            ReasonCode.LOW_COVERAGE,
            f"No called SNV position in {sample_name} was piled up: every site "
            "failed the pooled-coverage gate, or no spanning read belongs to a "
            "called cell.",
        )

    raw_snp = cell_counts.lazy()

    # Consensus support per called variant: Pos is 1-based amplicon-local in the
    # calls parquet, so -1 aligns it to the 0-based pileup coordinate; chrom is
    # the amplicon contig, which equals the per-cell `feature`. `callers` is the
    # sorted list of callers backing each allele `n_callers` is its length.
    callers_support = (
        pl.scan_parquet(calls_parquet_path)
        .filter(pl.col("type") == "SNV")
        .group_by(
            pl.col("chrom").alias(FEATURE),
            (pl.col("pos").cast(pl.Int32) - 1).alias(POS_LOCAL),
            pl.col("ref").alias(REF_OUT),
            pl.col("alt").alias(ALT_OUT),
        )
        .agg(
            pl.col("caller").unique().sort().alias(CALLERS),
            pl.col("caller").n_unique().cast(pl.Int32).alias(N_CALLERS),
        )
    )

    # Per-cell-support for variants: how many cells vote for each ALT, using the
    # same flag thresholds count_variant labels a per-cell call PASS with, grouped
    # on the same (feature, pos_local, ref, alt) as callers_support.
    cell_support = per_site_cell_support(
        raw_snp=raw_snp,
        flag_min_dp=args.flag_min_dp,
        flag_min_vaf=args.flag_min_vaf,
    )

    format_raw_snps_wide(
        raw_snp=raw_snp,
        callers_support=callers_support,
        cell_support=cell_support,
        flag_min_dp=args.flag_min_dp,
        flag_min_vaf=args.flag_min_vaf,
        flag_min_bq=args.flag_min_bq,
        flag_min_mq=args.flag_min_mq,
        genomic_map=genomic_map,
    ).sink_parquet("raw_variants.parquet")

    # Soft-filtered per-cell genotypes: one row per (barcode, cell-site), the
    # genotype and its GQ from the beta-binomial likelihood. Cells whose reads do
    # not separate the three genotypes are labelled lowGQ / lowqual.
    # Everything confident gets a real genotype, hom-ref included.
    all_snps = filter_snps(
        raw_snp=raw_snp,
        min_gq=args.min_gq,
        eps=args.error_eps,
        rho=args.error_rho,
    ).collect()

    # mod["snv"]: cells x the piled-up positions, X = non-reference VAF,
    # base composition and genotype in layers. Not every amplicon base: a
    # declared target the pileup never reached is in uns, not var.
    build_snv_anndata(
        raw_snp=raw_snp,
        genotypes=all_snps,
        calls=pl.scan_parquet(calls_parquet_path).filter(pl.col("type") == "SNV"),
        fasta_path=fasta_path,
    ).write_h5ad("variants_count.h5ad")

    # Carry consensus support onto each count row. callers_support is on
    # the same (feature, 0-based local pos, ref, alt) as it matches the raw parquet.
    # It is null for het/noise alleles no variant caller reported. The h5ad takes its own
    # site-level evidence straight from the catalog, so this join stays local to
    # the parquet.
    all_snps = all_snps.join(
        callers_support.collect(),
        on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT],
        how="left",
    )
    # Lift to genomic chrom/pos and generate schema order as before
    # (chrom/pos right after pos_local_0based; null for a plain FASTA).
    (
        add_genomic_coords(all_snps.lazy(), genomic_map)
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
        .sink_parquet("variants_count.parquet")
    )
else:
    # If no called positions for this sample still produce a schema-matched empty
    # raw_variants.parquet so the raw_variants_parquet output always has a file to read.
    empty_raw_variants().write_parquet("raw_variants.parquet")
    status.record_and_exit(
        ReasonCode.LOW_FEATURES,
        f"No SNVs were called for {sample_name}, so there is nothing to "
        "attribute per cell.",
    )
