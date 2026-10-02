process TABIX {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("*.tbi"), emit: tbi
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tabix --force --preset vcf ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: "\$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${vcf}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: "stub"
    END_VERSIONS
    """
}
