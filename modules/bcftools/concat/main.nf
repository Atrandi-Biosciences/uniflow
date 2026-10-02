process CONCAT_VCF {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcfs)

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Gather step for the scattered callers (FreeBayes, GATK HC). Intervals are
    // disjoint so concat with --allow-overlaps should be able to merge them back into a
    // single coordinate-sorted VCF identical to the whole-BAM call
    // Each per-interval VCF is force indexed first
    def prefix = "${meta.sample_name}.${caller}"
    """
    for vcf in ${vcfs}; do
        bcftools index --tbi --force \${vcf}
    done
    bcftools concat --allow-overlaps --output-type z --output ${prefix}.vcf.gz ${vcfs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.${caller}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
