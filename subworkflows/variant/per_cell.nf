// PER_CELL_VARIANTS — derive the merged per-cell variant count matrix from
// the consensus VCF in a single pileup pass.
//
// VCF_TO_PARQUET emits one long-format parquet per sample (rows: sample ×
// caller × variant; per-caller QUAL comes from INFO/VARQUAL[idx]). That
// parquet is the single catalog source. The per-cell pileup is
// caller-independent, so COUNT_VARIANT runs once per sample over the deduped
// union of every caller's consensus positions; which callers backed each
// variant is carried per row in the `callers` column rather than per-caller tracks.

include { VCF_TO_PARQUET         } from '../../modules/vcf_to_parquet/main'
include { MERGE_VARIANTS_PERCELL } from '../../modules/merge_variants_percell/main'
include { COUNT_VARIANT ; COUNT_VARIANT_INDEL ; DETECTION_POWER } from '../../modules/dna_processing.nf'


workflow PER_CELL_VARIANTS {
    take:
    consensus_vcf
    bam_with_index
    sample_barcode_mapping
    filtered_amplicon_reads
    amplicon_fasta

    main:
    VCF_TO_PARQUET(consensus_vcf)

    // One COUNT_VARIANT task per sample. count_variant*.py piles up the
    // deduped union of all positions across every caller in a
    // single pass and save this information in the `callers` column.
    VCF_TO_PARQUET.out.parquet
        .combine(bam_with_index, by: 0)
        .combine(sample_barcode_mapping, by: 0)
        .combine(filtered_amplicon_reads, by: 0)
        .combine(amplicon_fasta)
        .map { meta, parquet, bam, bai, barcode, filtered_reads, fasta ->
            tuple(meta, bam, barcode, fasta, parquet, bai, filtered_reads)
        }
        .set { count_variant_input }

    COUNT_VARIANT(count_variant_input)

    COUNT_VARIANT_INDEL(count_variant_input)

    // SNV only: the indel table borrows the SNV eps
    // An indel power number would restate this one with a known-wrong error
    // rate. The counts parquet is an optional emit, so a sample with no SNV
    // genotypes simply produces no report rather than an empty one.
    DETECTION_POWER(COUNT_VARIANT.out.variants_count_parquet)

    // The published per-cell tables: each pileup parquet merged with its
    // genotype sibling.
    COUNT_VARIANT.out.raw_variants_parquet
        .join(COUNT_VARIANT.out.variants_count_parquet, by: 0, remainder: true)
        .join(COUNT_VARIANT_INDEL.out.raw_variants_indel_parquet, by: 0)
        .join(COUNT_VARIANT_INDEL.out.variants_indel_count_parquet, by: 0, remainder: true)
        .map { meta, snv_raw, snv_counts, indel_raw, indel_counts ->
            tuple(meta, snv_raw, snv_counts ?: [], indel_raw, indel_counts ?: [])
        }
        .set { merge_input }

    MERGE_VARIANTS_PERCELL(merge_input)

    emit:
    variants_catalog_parquet = VCF_TO_PARQUET.out.parquet
    h5ad                     = COUNT_VARIANT.out.variants_count_h5ad
    counts_parquet           = COUNT_VARIANT.out.variants_count_parquet
    raw_parquet              = COUNT_VARIANT.out.raw_variants_parquet
    indel_h5ad               = COUNT_VARIANT_INDEL.out.variants_indel_count_h5ad
    indel_counts_parquet     = COUNT_VARIANT_INDEL.out.variants_indel_count_parquet
    indel_raw_parquet        = COUNT_VARIANT_INDEL.out.raw_variants_indel_parquet
    power_report             = DETECTION_POWER.out.report
    power_table              = DETECTION_POWER.out.table
    snv_parquet              = MERGE_VARIANTS_PERCELL.out.snv_parquet
    indel_parquet            = MERGE_VARIANTS_PERCELL.out.indel_parquet
}
