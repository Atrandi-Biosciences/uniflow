process AGGREGATE_METRICS {
    label 'small_job'
    tag "experiment_metrics"
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    debug false

    input:
    // stageAs deals with naming conflict
    tuple val(id), path(metrics_csv_files, stageAs: "?/*")

    output:
    tuple val(id), path("metrics.csv"), emit: metrics_csv

    script:
    // Build a safe, space-separated list of quoted paths
    """
        aggregate_metrics.py ${id} ${metrics_csv_files}
    """

    stub:
    """
    touch metrics.csv
    """
}
