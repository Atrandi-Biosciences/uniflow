process BCFTOOLS_FILTER {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.filtered.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Uniform discovery threshold on the per-caller normalised VCF (biallelic after
    // BCFTOOLS_NORM, so FORMAT/AD = [REF, ALT]).
    def expr = task.ext.args ?: 'FMT/AD[0:1] < 1'
    def prefix = "${meta.sample_name}.${caller}.filtered"
    """
    bcftools view -e '${expr}' -Oz -o ${prefix}.vcf.gz ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.${caller}.filtered"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
