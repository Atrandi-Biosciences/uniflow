include { ALIGN                 } from './align.nf'
include { FILTER                } from './filter.nf'
include { GET_SATURATION_CURVES } from '../../modules/rna_processing.nf'

workflow GENE_EXPRESSION {
    take:
    rna_libraries
    star_index_folder_ch
    whitelist_ch
    sample_barcode_mappings
    trimming_length
    internal_flag

    main:
    def modality_name = "gene_expression"

    ALIGN(
        rna_libraries,
        star_index_folder_ch,
        whitelist_ch,
        trimming_length,
        modality_name,
    )

    FILTER(
        ALIGN.out.sorted_bam,
        sample_barcode_mappings,
    )


    FILTER.out.to_saturation
        | combine(internal_flag)
        | GET_SATURATION_CURVES

    emit:
    metrics                       = channel.empty().mix(ALIGN.out.metrics).mix(GET_SATURATION_CURVES.out.saturation_csv).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    internal_metrics              = channel.empty().groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    qc                            = channel.empty().mix(FILTER.out.qc_plots).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    dev                           = channel.empty().mix(FILTER.out.merged_reads_parquet).mix(FILTER.out.raw_counts_df_parquet).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    out                           = channel.empty().mix(ALIGN.out.sorted_bam).mix(FILTER.out.bam_index).groupTuple(by: 0)
    metrics_parquet               = channel.empty().mix(FILTER.out.metrics_parquet, GET_SATURATION_CURVES.out.metrics_parquet)
    raw_gene_counts_h5ad          = FILTER.out.raw_gene_counts_h5ad
    cross_processing_top_cells    = FILTER.out.top_cells_parquet
    cross_processing_merged_reads = FILTER.out.merged_reads_parquet
}
