process GATK4_HAPLOTYPECALLER {
    tag "${meta.sample_name}.${interval_id}"
    label 'variant_calling_scatter'
    container 'community.wave.seqera.io/library/bcftools_gatk4_tabix:ddb88fdddba20b24'
    debug false

    input:
    tuple val(meta), path(bam), path(bai), val(interval_id), val(interval_line)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), val('gatkhc'), path("${meta.sample_name}.gatkhc.${interval_id}.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Scatter restricts this to a single amplicon interval (--intervals interval.bed).
    // GATK args (incl. --max-reads-per-alignment-start = the per-position DP cap)
    // come from conf/variant.config
    def args   = task.ext.args ?: ''
    def prefix = "${meta.sample_name}.gatkhc.${interval_id}"
    def avail_mem_mb = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    printf '%s\\n' "${interval_line}" > interval.bed
    gatk --java-options "-Xmx${avail_mem_mb}M -XX:+CrashOnOutOfMemoryError" HaplotypeCaller \\
        --reference ${fasta} \\
        --input ${bam} \\
        --output ${prefix}.vcf.gz \\
        --intervals interval.bed \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.gatkhc.${interval_id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "stub"
    END_VERSIONS
    """
}
