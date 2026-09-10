process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_gatk4_tabix:ddb88fdddba20b24'
    debug false

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk --java-options "-Xmx6g" CreateSequenceDictionary \\
        --REFERENCE ${fasta} \\
        --URI ${fasta} \\
        --TMP_DIR . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "stub"
    END_VERSIONS
    """
}
