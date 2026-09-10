include { EXTRACT_GENE_EXPRESSION_BAM ; FILTER_GENE_EXPRESSION_CELLS } from '../../modules/rna_processing.nf'
include { MERGE_READS                  } from '../../modules/shared.nf'

workflow FILTER {
    take:
    sorted_bam
    sample_barcode_mappings

    main:
    sorted_bam
        | EXTRACT_GENE_EXPRESSION_BAM

    EXTRACT_GENE_EXPRESSION_BAM.out.bam_df_parquet
        | combine(sample_barcode_mappings, by: [0])
        | MERGE_READS
        | FILTER_GENE_EXPRESSION_CELLS

    emit:
    to_saturation         = FILTER_GENE_EXPRESSION_CELLS.out.to_saturation
    qc_plots              = FILTER_GENE_EXPRESSION_CELLS.out.qc_plots
    raw_counts_df_parquet = FILTER_GENE_EXPRESSION_CELLS.out.raw_counts_df_parquet
    top_cells_parquet     = FILTER_GENE_EXPRESSION_CELLS.out.top_cells_parquet
    metrics_parquet       = FILTER_GENE_EXPRESSION_CELLS.out.metrics_parquet
    merged_reads_parquet  = MERGE_READS.out.merged_reads_parquet
    raw_gene_counts_h5ad  = FILTER_GENE_EXPRESSION_CELLS.out.raw_gene_counts_h5ad
    bam_index             = EXTRACT_GENE_EXPRESSION_BAM.out.index_bai
}
