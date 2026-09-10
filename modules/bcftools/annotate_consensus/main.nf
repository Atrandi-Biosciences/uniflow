process BCFTOOLS_ANNOTATE_CONSENSUS {
    tag "${meta.sample_name}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.sample_name}.consensus.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Strip every INFO carried over by bcftools merge, then synthesize the
    // consensus surface from the per-caller FORMAT columns. Each call to
    // +fill-tags adds one INFO summary using bcftools' built-in expressions.
    // The closing annotate restricts the FORMAT block to VAF / AD / DP / GT
    // so downstream tabular conversion (CSV / parquet / H5) sees a fixed
    // schema.
    //
    // FORMAT/VAF is computed from FORMAT/AD (REF, ALT) by +fill-tags's
    // built-in VAF rule. INFO/VAF_{min,max} expose the across-caller spread
    // (no median — bcftools' median over an even count returns a float that
    // doesn't fit a Number=1,Type=Integer field cleanly for the DP/AD
    // counterparts, so we report min/median/max for counts and min/max for
    // the fractional VAF).
    def prefix = "${meta.sample_name}.consensus"
    """
    bcftools annotate --remove 'INFO' ${vcf} \\
      | bcftools +fill-tags -- -t FORMAT/VAF \\
      | bcftools +fill-tags -- -t 'INFO/VAF_min:1=min(FORMAT/VAF)' \\
      | bcftools +fill-tags -- -t 'INFO/VAF_max:1=max(FORMAT/VAF)' \\
      | bcftools +fill-tags -- -t 'INFO/NCALLERS:1=N_PASS(FORMAT/DP>0)' \\
      | bcftools +fill-tags -- -t 'INFO/DP_median:1=int(median(FORMAT/DP))' \\
      | bcftools +fill-tags -- -t 'INFO/DP_max:1=int(max(FORMAT/DP))' \\
      | bcftools +fill-tags -- -t 'INFO/DP_min:1=int(min(FORMAT/DP))' \\
      | bcftools +fill-tags -- -t 'INFO/AD_REF_median:1=int(median(FORMAT/AD[:0,]))' \\
      | bcftools +fill-tags -- -t 'INFO/AD_REF_max:1=max(FORMAT/AD[:0,])' \\
      | bcftools +fill-tags -- -t 'INFO/AD_REF_min:1=min(FORMAT/AD[:0,])' \\
      | bcftools +fill-tags -- -t 'INFO/AD_ALT_median:1=int(median(FORMAT/AD[:1,]))' \\
      | bcftools +fill-tags -- -t 'INFO/AD_ALT_max:1=max(FORMAT/AD[:1,])' \\
      | bcftools +fill-tags -- -t 'INFO/AD_ALT_min:1=min(FORMAT/AD[:1,])' \\
      | bcftools annotate --remove '^INFO,^FORMAT/VAF,^FORMAT/AD,^FORMAT/DP,^FORMAT/GT,^FORMAT/FILTER,^FORMAT/VARQUAL' \\
        -Oz -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.sample_name}.consensus"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
