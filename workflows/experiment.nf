include { AGGREGATE_METRICS          } from '../modules/aggregation.nf'
include { GENERATE_MULTIQC_REPORT ; PREPARE_MULTIQC_REPORT; BUNDLE_QC_REPORT; } from '../modules/report.nf'


workflow EXPERIMENT {
    take:
    metrics_csv_ch
    metrics_reference_ch
    ch_multiqc_files
    multiqc_static_config_ch
    multiqc_css_ch
    internal_flag
    qc_merged_channels
    params_file

    main:
    metrics_csv_ch
        | map { _id, metrics_csv ->
            [params.experiment_id, metrics_csv]
        }
        | groupTuple(by: 0)
        | AGGREGATE_METRICS

    PREPARE_MULTIQC_REPORT(
        metrics_reference_ch,
        AGGREGATE_METRICS.out.metrics_csv,
        "aggregated_multiqc_config.yaml",
        "aggregated",
        internal_flag
    )

    ch_multiqc_files = ch_multiqc_files
        .mix(PREPARE_MULTIQC_REPORT.out.custom_metrics.flatten())
        .collect()
        .map { files ->
            ["experiment", files]
        }
    GENERATE_MULTIQC_REPORT(
        ch_multiqc_files,
        "experiment",
        PREPARE_MULTIQC_REPORT.out.multiqc_config,
        multiqc_static_config_ch,
        multiqc_css_ch,
    )

    qc_tarball_ch = GENERATE_MULTIQC_REPORT.out.report.combine(AGGREGATE_METRICS.out.metrics_csv.map{_id, metrics_csv -> [metrics_csv]}).combine(params_file)
    qc_flat = qc_merged_channels
        .flatMap { meta, modality, plots ->
            plots.flatten().collect { plot ->
                tuple(meta.sample_name, modality, plot)
            }
        }

    qc_bundle_input = qc_flat
        .map { sample_name, modality, plot ->
            tuple(
                sample_name,
                modality,
                plot.name,
                plot
            )
        }
        .collect(flat: false)
    BUNDLE_QC_REPORT(qc_tarball_ch, qc_bundle_input)

    emit:
    report                 = GENERATE_MULTIQC_REPORT.out.report
    aggregated_metrics_csv = AGGREGATE_METRICS.out.metrics_csv
    qc_tarball             = BUNDLE_QC_REPORT.out.tarball
}
