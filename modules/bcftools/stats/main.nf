process BCFTOOLS_STATS {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.bcftools_stats.txt"), emit: stats
    tuple val(meta), path("${meta.sample_name}.${caller}.bcftools_stats.txt"), emit: out_stats
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools stats ${vcf} > ${meta.sample_name}.${caller}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_name}.${caller}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
