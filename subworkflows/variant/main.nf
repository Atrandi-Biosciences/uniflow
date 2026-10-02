include { VARIANT_CALLING     } from './calling.nf'
include { VARIANT_HARMONIZE   } from './harmonize.nf'
include { VARIANT_MERGING     } from './merging.nf'
include { VARIANT_COMPARE     } from './compare.nf'
include { PER_CELL_VARIANTS   } from './per_cell.nf'
include { VARIANTS_TO_METRICS } from '../../modules/variants_to_metrics/main'
include { QC_VARIANTS         } from '../../modules/qc/qc_variants/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_CONSENSUS } from '../../modules/bcftools/stats/main'


workflow VARIANT {
    take:
    variant_calling_entry   // shared variant BAM -> bulk discovery + per-cell layer
    amplicon_fasta_file_ch
    sample_barcode_mappings
    filtered_amplicon_reads

    main:

    def modality_name = "variant"

    // Bulk discovery and the per-cell layer both read the same variant BAM; the
    // bulk callers cap per-position DP via their own CLI args (bulk_max_dp).
    VARIANT_CALLING(variant_calling_entry, amplicon_fasta_file_ch)

    VARIANT_HARMONIZE(VARIANT_CALLING.out.vcf)
    VARIANT_MERGING(VARIANT_HARMONIZE.out.harmonized)
    
    VARIANT_COMPARE(
        VARIANT_HARMONIZE.out.harmonized,
        VARIANT_CALLING.out.fasta,
        VARIANT_CALLING.out.fai,
    )

    PER_CELL_VARIANTS(
        VARIANT_MERGING.out.consensus_vcf,
        variant_calling_entry,
        sample_barcode_mappings,
        filtered_amplicon_reads,
        amplicon_fasta_file_ch,
    )

    // Per-caller metrics read the VCF-derived parquet; the per-site cell-support
    // metrics and the per-site per-cell-vs-abundance plot read the per-cell
    // raw_variants.parquet. raw_parquet is always emitted, so the join stays
    // 1:1 per sample and drives both the metrics and QC steps.
    PER_CELL_VARIANTS.out.variants_catalog_parquet
        | join(PER_CELL_VARIANTS.out.raw_parquet, by: 0)
        | set { variants_qc_inputs }
    VARIANTS_TO_METRICS(variants_qc_inputs)

    // QC_VARIANTS additionally renders the indel per-site per-cell-vs-abundance plot,
    // so it also takes the per-cell raw_variants_indel.parquet (4-tuple input).
    variants_qc_inputs
        | join(PER_CELL_VARIANTS.out.indel_raw_parquet, by: 0)
        | set { qc_variants_inputs }
    QC_VARIANTS(qc_variants_inputs)

    VARIANT_MERGING.out.consensus_vcf
        | map { meta, vcf -> tuple(meta, 'consensus', vcf) }
        | BCFTOOLS_STATS_CONSENSUS

    emit:
    metrics          = channel.empty()
    internal_metrics = channel.empty().mix(VARIANT_CALLING.out.out_stats).mix(BCFTOOLS_STATS_CONSENSUS.out.out_stats).mix(PER_CELL_VARIANTS.out.power_report).mix(PER_CELL_VARIANTS.out.power_table).map { meta, files ->
        return [meta, modality_name, files]
    }
    qc              = channel.empty().mix(QC_VARIANTS.out.png).mix(QC_VARIANTS.out.indel_png).mix(QC_VARIANTS.out.authenticity_png).mix(QC_VARIANTS.out.authenticity_indel_png).map { meta, files ->
        return [meta, modality_name, files]
    }
    dev             = channel.empty().mix(PER_CELL_VARIANTS.out.snv_parquet).mix(PER_CELL_VARIANTS.out.indel_parquet).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }

    out             = channel.empty().mix(VARIANT_MERGING.out.consensus_vcf).mix(VARIANT_MERGING.out.consensus_tbi)
    metrics_parquet = channel.empty().mix(VARIANTS_TO_METRICS.out.metrics_parquet).mix(QC_VARIANTS.out.metrics_parquet)
    h5ad            = channel.empty().mix(PER_CELL_VARIANTS.out.h5ad).mix(PER_CELL_VARIANTS.out.indel_h5ad)
}
