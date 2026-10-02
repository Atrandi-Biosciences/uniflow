// Read downsampling + per-cell coverage QC (COMB-555).
//
// Modules:
//   - PER_AMPLICON_COVERAGE_QC : pseudobulk per-amplicon depth (samtools coverage)
//   - PER_CELL_COVERAGE_QC     : per-(cell, amplicon) coverage tables + rank-knee plots
//   - SUBSAMPLE_KEEPLIST       : CB-aware per-(cell, amplicon) cap -> read-name keep-list + stats
//   - SUBSAMPLE_QC             : per-sample depth histogram with the cap marked, off those stats
//   - SUBSAMPLE_APPLY          : subset the aligned BAM to the keep-list (samtools view -N)

// pseudobulk per-amplicon coverage. `samtools coverage` gives one row per
// contig (contig == amplicon) with meandepth and coverage%. Diagnostic side
// output; per the design thread it does NOT and CANNOT reflect per-cell depth.
process PER_AMPLICON_COVERAGE_QC {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    debug false

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_amplicon_coverage.tsv"), emit: coverage
    tuple val(meta.sample_name), path("*_amplicon_coverage.tsv"), topic: per_sample_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools coverage ${bam} > ${meta.sample_name}_amplicon_coverage.tsv
    """

    stub:
    """
    touch ${meta.sample_name}_amplicon_coverage.tsv
    """
}

process PER_CELL_COVERAGE_QC {
    tag "${meta.sample_name}"
    label 'lazy_polars'
    container 'community.wave.seqera.io/library/uniflow:f4cea3af05b31ae3'
    debug false

    input:
    tuple val(meta), path(merged_reads), path(filtered_amplicon_reads), path(top_cells)

    output:
    tuple val(meta), path("*_per_cell_amplicon_depth.parquet"),        emit: dist_parquet
    tuple val(meta), path("*_per_cell_amplicon_depth_summary.txt"),    emit: dist_summary
    tuple val(meta), path("*.csv"),                                    emit: tables
    tuple val(meta), path("*_cap_impact_curve.png"),                   emit: cap_impact_plot
    tuple val(meta), path("*_rank_knee*.png"),                         emit: knee_plots

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    per_cell_coverage_qc.py ${merged_reads} ${filtered_amplicon_reads} ${top_cells} ${meta.sample_name} ${params.downsample_reads_per_cell_amplicon}
    """

    stub:
    """
    touch ${meta.sample_name}_per_cell_amplicon_depth.parquet
    touch ${meta.sample_name}_per_cell_amplicon_depth_summary.txt
    touch ${meta.sample_name}_cap_impact_sweep.csv
    touch ${meta.sample_name}_depth_histogram.csv
    touch ${meta.sample_name}_per_cell_amplicon_rank_knee_pooled.png
    touch ${meta.sample_name}_per_cell_amplicon_rank_knee_faceted.png
    touch ${meta.sample_name}_per_cell_total_rank_knee.png
    touch ${meta.sample_name}_cap_impact_curve.png
    """
}

// CB-aware per-(cell, amplicon). Computes a read-name keep-list
// (group_by(amplicon, CB) to sample(min(n, N))) and a per
// (cell, amplicon) stats table (n_reads_raw/n_reads_used/covered/capped). The
// keep-list is computed from full-depth inputs; the full-depth COUNT_AMPLICON
// matrix stays the canonical coverage signal.
process SUBSAMPLE_KEEPLIST {
    tag "${meta.sample_name}"
    label 'lazy_polars'
    container 'community.wave.seqera.io/library/pip_pandas_polars-bio_scanpy_scipy:32821fdc7ba9e09a'
    debug false

    input:
    tuple val(meta), path(merged_reads), path(filtered_amplicon_reads), path(top_cells)

    output:
    tuple val(meta), path("*_keep_reads.txt"),                              emit: keep_list
    tuple val(meta), path("*_per_cell_amplicon_subsample_stats.parquet"),   emit: stats

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    subsample_keeplist.py ${merged_reads} ${filtered_amplicon_reads} ${top_cells} ${meta.sample_name} ${params.downsample_reads_per_cell_amplicon} ${params.downsample_seed}
    """

    stub:
    """
    touch ${meta.sample_name}_keep_reads.txt
    touch ${meta.sample_name}_per_cell_amplicon_subsample_stats.parquet
    """
}

process SUBSAMPLE_QC {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/uniflow:f4cea3af05b31ae3'
    debug false

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_subsample_before_after.png"), emit: plot
    tuple val(meta.sample_name), path("*_subsample_before_after_mqc.png"), topic: per_sample_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_subsample_stats.py ${stats} ${meta.sample_name} ${params.downsample_reads_per_cell_amplicon}
    cp ${meta.sample_name}_subsample_before_after.png ${meta.sample_name}_subsample_before_after_mqc.png
    """

    stub:
    """
    touch ${meta.sample_name}_subsample_before_after.png
    touch ${meta.sample_name}_subsample_before_after_mqc.png
    """
}

// Subset the aligned BAM to the keep-list. `samtools view -N` is a
// lossless selection (mate/flags/tags preserved); no reconstruction or reheader.
process SUBSAMPLE_APPLY {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    debug false

    input:
    tuple val(meta), path(bam), path(keep_list)

    output:
    tuple val(meta), path("downsampled.bam"),     emit: bam
    tuple val(meta), path("downsampled.bam.bai"), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -b -N ${keep_list} -@ ${task.cpus} -o downsampled.bam ${bam}
    samtools index downsampled.bam
    """

    stub:
    """
    touch downsampled.bam downsampled.bam.bai
    """
}
