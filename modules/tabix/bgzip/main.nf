process BGZIP {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${vcf}.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bgzip --threads ${task.cpus} --force ${args} ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: "\$(echo \$(bgzip --version 2>&1) | sed 's/^.*bgzip //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${vcf}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: "stub"
    END_VERSIONS
    """
}
