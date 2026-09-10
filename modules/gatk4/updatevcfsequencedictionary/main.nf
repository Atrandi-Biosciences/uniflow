process GATK4_UPDATEVCFSEQUENCEDICTIONARY {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_gatk4_tabix:ddb88fdddba20b24'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf), path(bam), path(bai)

    output:
    tuple val(meta), val(caller), path("${meta.sample_name}.${caller}.reheader.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.sample_name}.${caller}.reheader"
    """
    gatk --java-options "-Xmx6g" IndexFeatureFile \\
        --input ${vcf} \\
        --tmp-dir .

    gatk --java-options "-Xmx6g" UpdateVCFSequenceDictionary \\
        --variant ${vcf} \\
        --source-dictionary ${bam} \\
        --output ${prefix}.vcf.gz \\
        --replace true \\
        --create-output-variant-index false \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.${caller}.reheader"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "stub"
    END_VERSIONS
    """
}
