process COUNT_AMPLICON {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/pip_polars_pydantic:927d9f10d3d028ec"
    debug false

    input:
    tuple val(meta), path(full_read_parquet), path(reference_fasta)

    output:

    tuple val(meta), path("filtered_amplicon_reads.parquet"), emit: filtered_amplicon_reads
    // Used for SNP filtering
    tuple val(meta), path("positions.parquet"), emit: positions_parquet
    tuple val(meta), path("positions.parquet"), path("amplicon_upset_data.parquet"), path("amplicon_counts.parquet"), emit: amplicon_qc
    // Both anchor tables are unpublished and stay in the work-dir for debugging only
    tuple val(meta), path("anchors.parquet"), emit: anchors_parquet
    tuple val(meta), path("*_anchors.csv"), emit: anchors_csv
    tuple val(meta), path("*_counts.parquet"), emit: filtered_amplicon_counts
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    script:

    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        count_amplicon.py ${full_read_parquet} ${reference_fasta} ${meta.sample_name} ${task.cpus}
        """

    stub:
    """
    touch filtered_amplicon_reads.parquet
    touch filtered_amplicon.h5ad
    touch positions.parquet
    touch amplicon_upset_data.parquet
    touch anchors.parquet
    touch ${meta.sample_name}_anchors.csv
    touch amplicon_counts.parquet
    touch metrics.parquet
    """
}


process QC_AMPLICON {
    label 'small_job'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug true

    input:
    tuple val(meta), path(positions), path(amplicon_upset_data), path(filtered_counts_parquet), path(top_cells_parquet)

    output:

    tuple val(meta), path("*.png"), emit: qc_plots_png
    tuple val(meta), path("*.csv"), emit: qc_metrics_csv
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    script:

    """
        qc_amplicon.py ${meta.sample_name} ${positions} ${amplicon_upset_data} ${filtered_counts_parquet} ${top_cells_parquet}
        """

    stub:
    """
    touch ${meta.sample_name}_qc.png
    touch ${meta.sample_name}_qc.csv
    touch metrics.parquet
    """
}

process COUNT_VARIANT {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(bam), path(barcode_parquet), path(reference_fasta), path(calls), path(index), path(filtered_amplicon_reads)

    output:
    tuple val(meta), path("variants_count.h5ad"), emit: variants_count_h5ad, optional: true
    tuple val(meta), path("variants_count.parquet"), emit: variants_count_parquet, optional: true
    tuple val(meta), path("raw_variants.parquet"), emit: raw_variants_parquet
    tuple val(meta), path("status.parquet"), topic: status, optional: true

    script:
    // Per-cell thresholds are named flags supplied via ext.args (conf/variant.config).
    def args = task.ext.args ?: ''
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        count_variant.py ${bam} ${reference_fasta} ${meta.sample_name} ${barcode_parquet} ${calls} ${filtered_amplicon_reads} ${args}
        """

    stub:
    """
    touch variants_count.h5ad
    touch variants_count.parquet
    touch raw_variants.parquet
    """
}

process COUNT_VARIANT_INDEL {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(bam), path(barcode_parquet), path(reference_fasta), path(calls), path(index), path(filtered_amplicon_reads)

    output:
    tuple val(meta), path("variants_indel_count.h5ad"), emit: variants_indel_count_h5ad, optional: true
    tuple val(meta), path("variants_indel_count.parquet"), emit: variants_indel_count_parquet, optional: true
    tuple val(meta), path("raw_variants_indel.parquet"), emit: raw_variants_indel_parquet
    tuple val(meta), path("status.parquet"), topic: status, optional: true

    script:
    def args = task.ext.args ?: ''
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        count_variant_indel.py ${bam} ${reference_fasta} ${meta.sample_name} ${barcode_parquet} ${calls} ${filtered_amplicon_reads} ${args}
        """

    stub:
    """
    touch variants_indel_count.h5ad
    touch variants_indel_count.parquet
    touch raw_variants_indel.parquet
    """
}

// Per-run detection power: Reads the genotype table COUNT_VARIANT just
// wrote plus the likelihood parameters it was called with, and reports what THIS
// run can detect: callability at the target GQ, the false-carrier budget, the
// ADO-corrected minimum detectable clone, the bulk ceiling in carrier cells, and
// which of those is binding. Always on, like the coverage QC; it is a QC artifact
// per run, not a static document.
process DETECTION_POWER {
    label 'small_job'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(variants_count_parquet)

    output:
    tuple val(meta), path("*_experimental_detection_power.txt"), emit: report
    tuple val(meta), path("*_experimental_detection_power.csv"), emit: table

    script:
    // Genotype-likelihood + bulk-floor values come from conf/variant.config, so
    // the report describes the model the run actually applied.
    def args = task.ext.args ?: ''
    """
    detection_power_report.py ${variants_count_parquet} ${meta.sample_name} ${args}
    """

    stub:
    """
    touch ${meta.sample_name}_experimental_detection_power.txt
    touch ${meta.sample_name}_experimental_detection_power.csv
    """
}

process EXTRACT_AMPLICON_BAM {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/pip_pandas_polars-bio_scanpy_scipy:32821fdc7ba9e09a"
    debug false

    input:
    tuple val(meta), path(bam)

    output:

    tuple val(meta), path("DNA_bam.parquet"), emit: dna_bam_parquet

    script:

    """
        extract_amplicon_bam.py ${bam} ${meta.sample_name} ${task.cpus}
        """

    stub:
    """
    touch DNA_bam.parquet
    """
}


process FILTER_AMPLICON_CELLS {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug false

    input:

    tuple val(meta), path(filtered_reads_parquet), path(full_read_parquet)

    output:
    tuple val(meta), path("*_top_cells.parquet"), emit: top_cells_parquet, optional: true
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet
    tuple val(meta), path("status.parquet"), topic: status, optional: true
    tuple val(meta.sample_name), path("*_knee_plot_mqc.png"), optional: true, topic: per_sample_report

    script:
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        filter_amplicon_cells.py ${filtered_reads_parquet} ${full_read_parquet} ${meta.sample_name} ${meta.force_cells}
        """

    stub:
    """
    touch ${meta.sample_name}_top_cells.parquet
    touch ${meta.sample_name}_qc.png
    touch metrics.parquet
    touch foo.h5ad
    """
}
