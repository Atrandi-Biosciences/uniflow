process VCF_TO_PARQUET {
    tag "${meta.sample_name}"
    label 'small_job'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.sample_name}.variants.parquet"), emit: parquet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf_to_parquet.py ${vcf} ${meta.sample_name} ${meta.sample_name}.variants.parquet
    """

    stub:
    """
    touch ${meta.sample_name}.variants.parquet
    """
}
