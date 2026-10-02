include { TRIM_READS } from '../../modules/shared.nf'
include { ALIGN_GENE_EXPRESSION } from '../../modules/rna_processing.nf'
workflow ALIGN {
    take:
        rna_libraries
        star_index_folder_ch
        whitelist_ch
        trimming_length
        modality_name
    main:

        TRIM_READS(
        rna_libraries,
        trimming_length,
        )
        ALIGN_GENE_EXPRESSION(
            TRIM_READS.out.trimmed_fastq,
            star_index_folder_ch,
            whitelist_ch,
            params.protocol,
            params.solofeature,
            modality_name,
        )
    emit:
        sorted_bam = ALIGN_GENE_EXPRESSION.out.bam_sorted
        metrics = channel.empty()
            .mix(ALIGN_GENE_EXPRESSION.out.log_final_txt)
            .mix(ALIGN_GENE_EXPRESSION.out.solo_metrics_csv)

}