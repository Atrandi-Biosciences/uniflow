include { GENE_EXPRESSION  } from '../subworkflows/gene_expression/main.nf'
include { CROSS_PROCESSING } from '../subworkflows/cross_processing/main.nf'

workflow RNA {
    take:
    rna_libraries
    star_index_folder_ch
    whitelist_ch
    sample_barcode_mappings
    trimming_length
    amplicon_fasta_file_ch
    run_gene_expression_ch
    run_cross_processing_ch

    main:
    // rna libraries if run_gene_expression_ch is true, else an empty channel
    rna_libraries.combine(run_gene_expression_ch)
        | map { meta, fastq1, fastq2, _run_gene_expression ->
            return [meta, fastq1, fastq2]
        }
        | set { gene_expression_entry }

    GENE_EXPRESSION(
        gene_expression_entry,
        star_index_folder_ch,
        whitelist_ch,
        sample_barcode_mappings,
        trimming_length,
    )

    rna_libraries.combine(run_cross_processing_ch)
        | map { meta, fastq1, fastq2, _run_cross_processing ->
            return [meta, fastq1, fastq2]
        }
        | set { cross_processing_entry }

    CROSS_PROCESSING(
        cross_processing_entry,
        amplicon_fasta_file_ch,
        GENE_EXPRESSION.out.cross_processing_merged_reads,
        GENE_EXPRESSION.out.cross_processing_top_cells,
    )

    emit:
    metrics          = channel.empty().mix(GENE_EXPRESSION.out.metrics).mix(CROSS_PROCESSING.out.metrics).groupTuple(by: [0, 1])
    internal_metrics = channel.empty().mix(GENE_EXPRESSION.out.internal_metrics).mix(CROSS_PROCESSING.out.internal_metrics).groupTuple(by: [0, 1])
    qc               = channel.empty().mix(GENE_EXPRESSION.out.qc).mix(CROSS_PROCESSING.out.qc).groupTuple(by: [0, 1])
    dev              = channel.empty().mix(GENE_EXPRESSION.out.dev).mix(CROSS_PROCESSING.out.dev).groupTuple(by: [0, 1])
    out              = channel.empty().mix(GENE_EXPRESSION.out.out).mix(CROSS_PROCESSING.out.out).groupTuple(by: 0)
    metrics_parquet  = channel.empty().mix(GENE_EXPRESSION.out.metrics_parquet).mix(CROSS_PROCESSING.out.metrics_parquet)
    counts_h5ad      = channel.empty().mix(GENE_EXPRESSION.out.raw_gene_counts_h5ad)
}
