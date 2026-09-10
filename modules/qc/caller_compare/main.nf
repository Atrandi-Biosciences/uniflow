process CALLER_COMPARE_QC {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(callers), path(vcfs)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.sample_name}_caller_compare/"), emit: report
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Stage with caller-only filenames so compare_caller_runs.sh derives
    // labels = caller (it strips .vcf.gz / .norm). VCFs come in caller order.
    def n = callers.size()
    def stage = (0..<n).collect { i -> "cp ${vcfs[i]} ${callers[i]}.vcf.gz" }.join('\n    ')
    def vcf_args = callers.collect { "${it}.vcf.gz" }.join(' ')
    def outdir = "${meta.sample_name}_caller_compare"
    """
    ${stage}

    compare_caller_runs.sh \\
        --ref ${fasta} \\
        --out ${outdir} \\
        ${vcf_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def outdir = "${meta.sample_name}_caller_compare"
    """
    mkdir -p ${outdir}
    touch ${outdir}/summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
