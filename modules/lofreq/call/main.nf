process LOFREQ_CALL {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4'
    debug false

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val('lofreq'), path("${meta.sample_name}.lofreq.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.sample_name}.lofreq"
    """
    lofreq call-parallel \\
        --pp-threads ${task.cpus} \\
        ${args} \\
        -f ${fasta} \\
        -o ${prefix}.vcf.gz \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: "\$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.lofreq"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: "stub"
    END_VERSIONS
    """
}
