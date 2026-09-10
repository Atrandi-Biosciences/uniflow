process SAMTOOLS_FAIDX {
    label 'small_job'
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    debug false

    input:
    path fasta

    output:
    path fasta, emit: fasta
    path "${fasta}.fai", emit: fai
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}
