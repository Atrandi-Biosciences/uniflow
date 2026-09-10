// SAMTOOLS_SORT — coordinate-sort a BAM and write its index, emitting (meta, bam, bai).
//
// Single-job module, currently unused. Kept as the sort half of the deferred samtools split
// out of the ALIGN_AMPLICON aligner (see memory samtools-module-split): the aligner's inline
// `samtools sort` / `flagstats` should move into their own one-job modules; this is ready to
// wire in then.

process SAMTOOLS_SORT {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    debug false

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample_name}.sorted.bam"), path("${meta.sample_name}.sorted.bam.bai"), emit: bam_with_index
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.sample_name}.sorted"
    """
    samtools sort -@ ${task.cpus - 1} ${args} \\
        --write-index -o ${prefix}.bam##idx##${prefix}.bam.bai ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_name}.sorted.bam
    touch ${meta.sample_name}.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}
