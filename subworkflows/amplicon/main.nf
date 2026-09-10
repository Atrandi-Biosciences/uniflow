include { COUNT_AMPLICON ; QC_AMPLICON ; EXTRACT_AMPLICON_BAM ; FILTER_AMPLICON_CELLS } from '../../modules/dna_processing.nf'
include { ALIGN_AMPLICON           } from '../../modules/minimap2/main'
include { MERGE_READS              } from '../../modules/shared.nf'
include { PER_AMPLICON_COVERAGE_QC ; PER_CELL_COVERAGE_QC ; SUBSAMPLE_KEEPLIST ; SUBSAMPLE_QC ; SUBSAMPLE_APPLY } from '../../modules/downsample/main.nf'


workflow AMPLICON {
    take:
    amplicon_entry
    amplicon_fasta_file_ch
    sample_barcode_mappings

    main:
    def modality_name = "amplicon"

    ALIGN_AMPLICON(
        amplicon_entry.map { meta, read1, _read2 -> [meta, read1] },
        amplicon_fasta_file_ch,
    )

    // Downsampling (per_cell_per_amplicon_downsample, on by default):
    // false => pass-through (full-depth variant-layer BAM).
    // true  => CB-aware per-(cell, amplicon) cap; the cap is computed from
    //          full-depth counts/cells, so counts and cells stay at full depth
    //          here and only the variant-layer BAM is capped.
    def amplicon_bam = ALIGN_AMPLICON.out.bam
    def amplicon_index = ALIGN_AMPLICON.out.index

    if (params.per_cell_per_amplicon_downsample) {
        def cap = params.downsample_reads_per_cell_amplicon
        if (cap == null || !(cap instanceof Number) || cap <= 0) {
            error("per_cell_per_amplicon_downsample requires downsample_reads_per_cell_amplicon > 0 (got: ${cap}).")
        }
    }

    amplicon_bam
        | EXTRACT_AMPLICON_BAM
        | join(sample_barcode_mappings, by: [0])
        | MERGE_READS

    MERGE_READS.out.merged_reads_parquet
        | combine(amplicon_fasta_file_ch)
        | COUNT_AMPLICON

    COUNT_AMPLICON.out.filtered_amplicon_counts
        | combine(MERGE_READS.out.merged_reads_parquet, by: 0)
        | FILTER_AMPLICON_CELLS

    COUNT_AMPLICON.out.amplicon_qc
        | join(FILTER_AMPLICON_CELLS.out.top_cells_parquet, by: 0)
        | QC_AMPLICON

    amplicon_bam
        | combine(amplicon_index, by: 0)
        | set { aligned_bam_with_index }

    // CB-aware per-(cell, amplicon) cap. The keep-list is computed from the
    // full-depth read->amplicon join (COUNT_AMPLICON), read->CB (MERGE_READS)
    // and called cells (FILTER_AMPLICON_CELLS); `samtools view -N` then subsets
    // the full aligned BAM losslessly. Counts/cells above stay full-depth; only
    // the variant-layer BAM (discovery + per-cell) is capped.
    def variant_bam_bai_ch = aligned_bam_with_index
    def subsample_qc_plot = channel.empty()
    if (params.per_cell_per_amplicon_downsample) {
        MERGE_READS.out.merged_reads_parquet
            | combine(COUNT_AMPLICON.out.filtered_amplicon_reads, by: 0)
            | combine(FILTER_AMPLICON_CELLS.out.top_cells_parquet, by: 0)
            | SUBSAMPLE_KEEPLIST

        amplicon_bam
            | combine(SUBSAMPLE_KEEPLIST.out.keep_list, by: 0)
            | SUBSAMPLE_APPLY

        variant_bam_bai_ch = SUBSAMPLE_APPLY.out.bam
            | combine(SUBSAMPLE_APPLY.out.index, by: 0)

        // The stats parquet is not published; SUBSAMPLE_QC renders it instead.
        SUBSAMPLE_KEEPLIST.out.stats | SUBSAMPLE_QC
        subsample_qc_plot = SUBSAMPLE_QC.out.plot
    }

    PER_AMPLICON_COVERAGE_QC(amplicon_bam)

    MERGE_READS.out.merged_reads_parquet
        | combine(COUNT_AMPLICON.out.filtered_amplicon_reads, by: 0)
        | combine(FILTER_AMPLICON_CELLS.out.top_cells_parquet, by: 0)
        | PER_CELL_COVERAGE_QC

    def coverage_qc_metrics = PER_AMPLICON_COVERAGE_QC.out.coverage
        .mix(PER_CELL_COVERAGE_QC.out.dist_summary)
        .mix(PER_CELL_COVERAGE_QC.out.tables)
    def coverage_qc_plots = PER_CELL_COVERAGE_QC.out.plots.mix(subsample_qc_plot)

    emit:
    metrics                         = channel.empty().mix(ALIGN_AMPLICON.out.mapping_rates).mix(coverage_qc_metrics).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    internal_metrics                = channel.empty().mix(QC_AMPLICON.out.qc_metrics_csv).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    qc                              = channel.empty().mix(QC_AMPLICON.out.qc_plots_png).mix(FILTER_AMPLICON_CELLS.out.qc_plots).mix(coverage_qc_plots).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    dev                             = channel.empty().mix(COUNT_AMPLICON.out.filtered_amplicon_reads).mix(COUNT_AMPLICON.out.filtered_amplicon_counts).mix(MERGE_READS.out.merged_reads_parquet).mix(PER_CELL_COVERAGE_QC.out.dist_parquet).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    out                             = channel.empty().mix(amplicon_bam).mix(amplicon_index).groupTuple(by: 0)
    h5ad                            = FILTER_AMPLICON_CELLS.out.h5ad
    metrics_parquet                 = channel.empty().mix(COUNT_AMPLICON.out.metrics_parquet).mix(QC_AMPLICON.out.metrics_parquet).mix(FILTER_AMPLICON_CELLS.out.metrics_parquet)
    variant_bam_bai                 = variant_bam_bai_ch
    variant_filtered_amplicon_reads = COUNT_AMPLICON.out.filtered_amplicon_reads
}
