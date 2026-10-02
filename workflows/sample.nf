include { MERGE_H5AD ; MERGE_METRICS as MERGE_SAMPLE_METRICS ; MERGE_STATUS } from '../modules/shared.nf'
include { DOWNSTREAM_ANALYSIS                   } from '../modules/downstream_analysis.nf'
include { GENERATE_MULTIQC_REPORT ; PREPARE_MULTIQC_REPORT } from '../modules/report.nf'

workflow SAMPLE {
    take:
    h5ad_files
    ch_sample_metrics_files
    counts_format
    metrics_reference_ch
    multiqc_static_config_ch
    multiqc_css_ch
    internal_flag

    main:
    h5ad_files
        | map { meta, h5ad ->
            [meta - meta.subMap("library_type", "library_id"), h5ad]
        }
        | groupTuple(by: 0)
        | map { meta, h5ads ->
            [meta, file(params.input_csv), params.experiment_id, h5ads]
        }
        | MERGE_H5AD

    def counts_output = channel.empty()
    if (counts_format.contains("h5ad")) {
        counts_output = counts_output.mix(MERGE_H5AD.out.counts_h5ad)
    }
    if (counts_format.contains("zarr")) {
        counts_output = counts_output.mix(MERGE_H5AD.out.counts_zarr)
    }
    counts_output = counts_output
        | groupTuple(by: 0)

    ch_sample_metrics_files
        | map { meta, metrics ->
            [meta.sample_name, metrics]
        }
        | groupTuple(by: 0)
        | MERGE_SAMPLE_METRICS

    channel.topic("status")
        | map { meta, status ->
            [meta.sample_name, status]
        }
        | groupTuple(by: 0)
        | MERGE_STATUS
    PREPARE_MULTIQC_REPORT(
        metrics_reference_ch,
        MERGE_SAMPLE_METRICS.out.metrics_csv,
        "aggregated_multiqc_config.yaml",
        "single",
        internal_flag,
    )

    ch_multiqc_files = PREPARE_MULTIQC_REPORT.out.custom_metrics_with_id
    // Need to normalize the multiqc files into individual tuples for further processing
    ch_multiqc_files = ch_multiqc_files 
        | flatMap { sample_id, files ->
            files.collect { file ->
                tuple(sample_id, file)
            }
        }
    per_sample_report = channel.topic("per_sample_report")
    ch_multiqc_files = ch_multiqc_files.mix(per_sample_report).groupTuple(by: 0)
    GENERATE_MULTIQC_REPORT(
        ch_multiqc_files,
        "report",
        PREPARE_MULTIQC_REPORT.out.multiqc_config,
        multiqc_static_config_ch,
        multiqc_css_ch,
    )

    emit:
    counts             = counts_output
    per_sample_metrics = MERGE_SAMPLE_METRICS.out.metrics_csv
    per_sample_status  = MERGE_STATUS.out.status_parquet
    per_sample_report  = GENERATE_MULTIQC_REPORT.out.report_with_id
}
