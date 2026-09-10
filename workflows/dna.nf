include { AMPLICON } from '../subworkflows/amplicon/main.nf'
include { VARIANT  } from '../subworkflows/variant/main.nf'

workflow DNA {
    take:
    dna_libraries
    amplicon_fasta_file_ch
    sample_barcode_mappings
    run_amplicon_ch
    run_variant_calling_ch

    main:

    // Amplicon modality start
    dna_libraries.combine(run_amplicon_ch)
        | map { meta, fastq1, fastq2, _run_amplicon ->
            return [meta, fastq1, fastq2]
        }
        | set { amplicon_entry }

    AMPLICON(
        amplicon_entry,
        amplicon_fasta_file_ch,
        sample_barcode_mappings,
    )
    // Amplicon modality end

    // Variant modality start
    AMPLICON.out.variant_bam_bai.combine(run_variant_calling_ch)
        | map { meta, bam, bai, _run_variant_calling ->
            return [meta, bam, bai]
        }
        | set { variant_calling_entry }

    VARIANT(
        variant_calling_entry,
        amplicon_fasta_file_ch,
        sample_barcode_mappings,
        AMPLICON.out.variant_filtered_amplicon_reads,
    )

    emit:
    metrics          = channel.empty().mix(AMPLICON.out.metrics).mix(VARIANT.out.metrics).groupTuple(by: [0, 1])
    internal_metrics = channel.empty().mix(AMPLICON.out.internal_metrics).mix(VARIANT.out.internal_metrics).groupTuple(by: [0, 1])
    qc               = channel.empty().mix(AMPLICON.out.qc).mix(VARIANT.out.qc).groupTuple(by: [0, 1])
    dev              = channel.empty().mix(AMPLICON.out.dev).mix(VARIANT.out.dev).groupTuple(by: [0, 1])
    out              = channel.empty().mix(AMPLICON.out.out).mix(VARIANT.out.out).groupTuple(by: 0)
    metrics_parquet  = channel.empty().mix(AMPLICON.out.metrics_parquet).mix(VARIANT.out.metrics_parquet)
    counts_h5ad      = channel.empty().mix(AMPLICON.out.h5ad).mix(VARIANT.out.h5ad)
}
