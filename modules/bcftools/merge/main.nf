process BCFTOOLS_MERGE {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `--info-rules VARQUAL:join` concatenates each input file's
    // INFO/VARQUAL into a per-record list
    args   = task.ext.args   ?: '--force-samples -m none --info-rules VARQUAL:join'
    prefix = task.ext.prefix ?: "${meta.sample_name}.merged"
    def vcf_list = vcfs.collect { it.toString() }.join(' ')
    """
    for v in ${vcf_list}; do
        bcftools index --force --tbi \$v
    done

    bcftools merge ${args} ${vcf_list} -Oz -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.sample_name}.merged"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
