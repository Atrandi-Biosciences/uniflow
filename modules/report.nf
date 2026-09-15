#!/usr/bin/env nextflow

process GENERATE_MULTIQC_REPORT {

    container "community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c"
    debug false

    input:
    tuple val(source_id), path("*", stageAs: "?/*")
    val output_name
    path aggregated_config_path
    path static_config, stageAs: "static_multiqc_config.yaml"
    path custom_css, stageAs: "uniflow_report.css"

    output:
    path "${output_name}.html", emit: report
    tuple val(source_id), path("${output_name}.html"), emit: report_with_id
    path "${output_name}_data", emit: data

    script:
    """
    multiqc -e general_stats . -n ${output_name}.html -c ${static_config} -c ${aggregated_config_path}
    """

    stub:
    """
    touch ${output_name}.html
    mkdir -p ${output_name}_data
    """
}

process PREPARE_MULTIQC_REPORT {
    container "community.wave.seqera.io/library/pip_pyyaml_multiqc_polars:d0a8456de2ce4dce"
    debug false

    input:
    path metrics_reference
    tuple val(id), path(aggregated_metrics)
    val output
    val scope
    val internal_flag

    output:
    path "${output}", emit: multiqc_config
    path "*.tsv", emit: custom_metrics
    tuple val(id), path("*.tsv"), emit: custom_metrics_with_id

    script:
    """
    prepare_multiqc_report.py ${metrics_reference} ${id} ${aggregated_metrics} ${scope} ${output} ${internal_flag}
    """

    stub:
    """
    touch ${output}
    touch ${scope}_custom_metrics.tsv
    """
}

process BUNDLE_QC_REPORT {
    input:
    tuple path(experiment_metrics_csv),
          path(aggregated_params),
          path(aggregated_multiqc_report)

    val qc_entries

    output:
    path "qc_report.tar.gz", emit: tarball

    script:
    def copyCommands = qc_entries.collect { sample_name, modality, filename, plot ->
        """
        mkdir -p plots/${sample_name}/${modality}
        cp -L "${plot}" "plots/${sample_name}/${modality}/${filename}"
        """
    }.join('\n')

    """
    ${copyCommands}

    tar -chvzf qc_report.tar.gz \
        ${experiment_metrics_csv} \
        ${aggregated_params} \
        ${aggregated_multiqc_report} \
        plots
    """

    stub:
    """
    touch qc_report.tar.gz
    """
}