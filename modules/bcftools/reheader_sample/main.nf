process BCFTOOLS_REHEADER_SAMPLE {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.reheader.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Single-sample input. Rename the sole sample column to the caller string
    // so downstream `bcftools merge` gives us caller-as-pseudo-sample layout.
    def prefix = "${meta.sample_name}.${caller}.reheader"
    // Strip INFO/TYPE, INFO/MQ, FORMAT/GQ here: callers declare these with
    // conflicting Number/Type in their headers, which makes bcftools merge
    // abort once a record's value doesn't fit the picked definition.
    """
    echo "${caller}" > samples.txt
    bcftools reheader --samples samples.txt ${vcf} \\
        | bcftools annotate -x INFO/TYPE,INFO/MQ,FORMAT/GQ \\
        | bcftools view -Oz -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.${caller}.reheader"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
