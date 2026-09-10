process BEDTOOLS_MAKEWINDOWS {
    label 'small_job'
    container 'community.wave.seqera.io/library/bedtools:2.31.1--db419eb70de48d65'
    debug false

    input:
    path fasta_fai

    output:
    path "${fasta_fai.simpleName}.bed", emit: bed
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bedtools makewindows -g ${fasta_fai} -n 1 > ${fasta_fai.simpleName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: "\$(bedtools --version | sed 's/^bedtools v//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta_fai.simpleName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: "stub"
    END_VERSIONS
    """
}
