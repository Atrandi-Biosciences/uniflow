process VARDICT {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'community.wave.seqera.io/library/bcftools_tabix_vardict-java:18655879b49d9e3b'
    debug false

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path amplicon_bed

    output:
    tuple val(meta), val('vardict'), path("${meta.sample_name}.vardict.vcf"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // VarDict args (floors + BED column indices) come from conf/variant.config
    def args   = task.ext.args  ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = "${meta.sample_name}.vardict"
    """
    vardict-java \\
        ${args} \\
        -th ${task.cpus} \\
        -G ${fasta} \\
        -b ${bam} \\
        ${amplicon_bed} \\
        | teststrandbias.R \\
        | var2vcf_valid.pl \\
            ${args2} \\
            -N ${meta.sample_name} \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: "\$(realpath \$(command -v vardict-java) | sed 's/.*java-//; s/-.*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.vardict"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: "stub"
    END_VERSIONS
    """
}
