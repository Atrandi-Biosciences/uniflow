include { AGGREGATE_METRICS          } from '../modules/aggregation.nf'
include { GENERATE_AGGREGATED_REPORT ; GENERATE_MULTIQC_CUSTOM } from '../modules/report.nf'


workflow EXPERIMENT {
    take:
    metrics_csv_ch
    metrics_reference_ch
    ch_multiqc_files
    multiqc_static_config_ch
    multiqc_css_ch

    main:
    metrics_csv_ch
        | map { _id, metrics_csv ->
            [params.experiment_id, metrics_csv]
        }
        | groupTuple(by: 0)
        | AGGREGATE_METRICS

    GENERATE_MULTIQC_CUSTOM(
        metrics_reference_ch,
        AGGREGATE_METRICS.out.metrics_csv,
        "aggregated_multiqc_config.yaml",
        "aggregated",
    )

    ch_multiqc_files = ch_multiqc_files
        .mix(GENERATE_MULTIQC_CUSTOM.out.custom_metrics.flatten())
        .collect()
        .map { files ->
            ["experiment", files]
        }
    GENERATE_AGGREGATED_REPORT(
        ch_multiqc_files,
        "experiment",
        GENERATE_MULTIQC_CUSTOM.out.multiqc_config,
        multiqc_static_config_ch,
        multiqc_css_ch,
    )

    emit:
    report                 = GENERATE_AGGREGATED_REPORT.out.report
    aggregated_metrics_csv = AGGREGATE_METRICS.out.metrics_csv
}
