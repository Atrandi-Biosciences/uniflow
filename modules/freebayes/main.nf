process FREEBAYES {
    tag "${meta.sample_name}.${interval_id}"
    label 'freebayes'
    container 'community.wave.seqera.io/library/bcftools_freebayes_tabix:0ab6fc35c74cce4b'
    debug false

    input:
    tuple val(meta), path(bam), path(bai), val(interval_id), val(interval_line)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val('freebayes'), path("${meta.sample_name}.freebayes.${interval_id}.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Scatter restricts this to a single amplicon interval (--targets interval.bed).
    // FreeBayes args (floors + the --limit-coverage DP cap) come from
    // conf/variant.config
    def args = task.ext.args ?: ''
    def prefix = "${meta.sample_name}.freebayes.${interval_id}"
    """
    printf '%s\\n' "${interval_line}" > interval.bed
    freebayes ${args} --targets interval.bed --fasta-reference ${fasta} ${bam} > ${prefix}.vcf
    bgzip --force ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: "\$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.freebayes.${interval_id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: "stub"
    END_VERSIONS
    """
}
