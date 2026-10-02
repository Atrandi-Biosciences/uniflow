process VARIANTS_TO_METRICS {
    tag "${meta.sample_name}"
    label 'small_job'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(parquet), path(raw_parquet)

    output:
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    variants_to_metrics.py ${parquet} ${meta.sample_name} metrics.parquet ${raw_parquet}
    """

    stub:
    """
    touch metrics.parquet
    """
}
