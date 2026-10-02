process BCFTOOLS_INFO_TO_FORMAT {
    tag "${meta.sample_name}.${caller}"
    label 'small_job'
    container 'community.wave.seqera.io/library/bcftools_tabix:4adcae3729d5f19f'
    debug false

    input:
    tuple val(meta), val(caller), path(vcf)

    output:
    tuple val(meta), val(caller), path("${prefix}.vcf.gz"), emit: vcf
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Lift INFO/FILTER and INFO/VARQUAL into the per-caller FORMAT block so
    // the sample-column-equals-caller scheme carries filter and quality per
    // caller through the downstream merge. `bcftools annotate -c FORMAT/X`
    // requires the new FORMAT fields to be declared in the header, so we
    // emit hdr.txt inline.
    prefix = task.ext.prefix ?: "${meta.sample_name}.${caller}.format_lifted"
    """
    cat <<-HDR > hdr.txt
    ##FORMAT=<ID=FILTER,Number=1,Type=String,Description="FILTER status">
    ##FORMAT=<ID=VARQUAL,Number=1,Type=Float,Description="Variant call quality">
    HDR

    bcftools query -f '%CHROM\\t%POS\\t%INFO/FILTER\\t%INFO/VARQUAL\\n' ${vcf} \\
        | bgzip -c > annot.txt.gz
    tabix -s1 -b2 -e2 annot.txt.gz

    bcftools annotate -a annot.txt.gz -h hdr.txt \\
        -c CHROM,POS,FORMAT/FILTER,FORMAT/VARQUAL \\
        ${vcf} -Oz -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.sample_name}.${caller}.format_lifted"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
