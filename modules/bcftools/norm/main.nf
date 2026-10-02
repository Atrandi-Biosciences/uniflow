process BCFTOOLS_NORM {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.norm.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Normalizes (split multiallelics + left-align). QUAL/DP filtering is not part of this step.
    // The only hard drop here is the spanning-deletion star allele, which carries no genotype of its own.
    def args   = task.ext.args ?: '--multiallelics -both --rm-dup all'
    def prefix = "${meta.sample_name}.${caller}.norm"
    """
    bcftools norm ${args} --fasta-ref ${fasta} ${vcf} \\
        | bcftools view -e 'ALT="*"' -Oz -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.${caller}.norm"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
