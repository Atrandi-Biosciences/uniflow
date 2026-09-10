process TEMPLATE_VCF {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    val(meta)

    output:
    tuple val(meta),
          path("${meta.sample_name}.template.vcf.gz"),
          path("${meta.sample_name}.template.vcf.gz.tbi"),
          emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Header-only VCF with a declared sample column. Used as the right-hand
    // side of `bcftools merge` to give a site-only VCF (e.g. lofreq) a sample
    // column so downstream `bcftools reheader --samples` can rename it.
    def prefix = "${meta.sample_name}.template"
    def sample = meta.sample_name
    """
    {
        echo '##fileformat=VCFv4.2'
        echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t%s\\n' '${sample}'
    } | bgzip -c > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.template"
    """
    touch ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
