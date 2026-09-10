process BCFTOOLS_CALL {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val('bcftools'), path("${meta.sample_name}.bcftools.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Args (the FORMAT/AD annotations the downstream AD floor needs + the
    // --max-depth per-position cap) come from conf/variant.config, the single source
    // of truth for the thresholds. The fallback keeps the required annotations only.
    def mpileup_args = task.ext.args  ?: '--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/SCR'
    def call_args    = task.ext.args2 ?: '--multiallelic-caller --variants-only --ploidy 2'
    def prefix = "${meta.sample_name}.bcftools"
    """
    bcftools mpileup ${mpileup_args} --fasta-ref ${fasta} ${bam} \\
        | bcftools call ${call_args} --output-type u \\
        | bcftools +fill-tags --output-type z --output ${prefix}.vcf.gz -- -t "FORMAT/VAF"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.bcftools"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
