#!/usr/bin/env nextflow

process FALCO {
    label 'small_job'
    container "community.wave.seqera.io/library/falco:1.2.5--bd70431b39ecac2f"
    debug false

    input:
    tuple val(library_id), val(library_type), path(read), val(validation_done)

    output:
    path ("*_fastqc_data.txt"), emit: txt
    path ("*_fastqc_report.html"), emit: html
    tuple val(library_id), val(library_type), path("*_fastqc_report.html"), emit: tagged_html
    path "versions.yml", topic: versions

    script:
    """
    falco ${read} \\
        --threads ${task.cpus} \\
        -data-filename ${read}_fastqc_data.txt \\
        -report-filename ${read}_fastqc_report.html \\
        -summary-filename ${read}_summary.txt

    for f in *_fastqc_data.txt; do
        [[ -e "\$f" ]] || continue
        sed -i '1s/^##Falco/##FastQC/' "\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: "\$(falco --version | sed -e "s/falco //g")"
    END_VERSIONS
    """

    stub:
    """
    touch ${read}_fastqc_data.txt
    touch ${read}_fastqc_report.html
    touch ${read}_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: "stub"
    END_VERSIONS
    """
}
