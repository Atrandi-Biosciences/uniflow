process QC_VARIANTS {
    tag "${meta.sample_name}"
    label 'small_job'
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug false

    input:
    tuple val(meta), path(parquet), path(raw_parquet), path(indel_raw_parquet)

    output:
    tuple val(meta), path("${meta.sample_name}_caller_upset.png"), emit: png
    tuple val(meta), path("${meta.sample_name}_caller_upset_indel.png"), emit: indel_png
    tuple val(meta), path("${meta.sample_name}_persite_percell_vs_abundance.png"), emit: authenticity_png
    tuple val(meta), path("${meta.sample_name}_persite_percell_vs_abundance_indel.png"), emit: authenticity_indel_png
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    qc_variants.py \\
        ${parquet} \\
        ${meta.sample_name} \\
        ${meta.sample_name}_caller_upset.png \\
        metrics.parquet \\
        ${raw_parquet} \\
        ${meta.sample_name}_persite_percell_vs_abundance.png \\
        ${indel_raw_parquet} \\
        ${meta.sample_name}_persite_percell_vs_abundance_indel.png \\
        ${meta.sample_name}_caller_upset_indel.png
    """

    stub:
    """
    touch ${meta.sample_name}_caller_upset.png
    touch ${meta.sample_name}_caller_upset_indel.png
    touch ${meta.sample_name}_persite_percell_vs_abundance.png
    touch ${meta.sample_name}_persite_percell_vs_abundance_indel.png
    touch metrics.parquet
    """
}
