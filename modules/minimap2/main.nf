// ALIGN_AMPLICON — single-mode minimap2 short-read aligner for the amplicon BAM.
//
// NOTE (deferred, see memory samtools-module-split): the inline `samtools sort` and
// `samtools flagstats` should later move into their own one-job modules.

process ALIGN_AMPLICON {
    tag "${meta.sample_name}"
    label 'cpu_low_mem'
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    debug false

    input:
    tuple val(meta), path(read1)
    path reference

    output:
    tuple val(meta), path("aligned.sorted.bam"),     emit: bam
    tuple val(meta), path("aligned.sorted.bam.bai"), emit: index
    tuple val(meta), path("*_mapping_rates.txt"),    emit: mapping_rates
    path "versions.yml", topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def rg = "@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:ILLUMINA"
    """
    minimap2 -ax sr --secondary=no --frag=no -R '${rg}' -t ${task.cpus} ${args} ${reference} ${read1} \\
      | samtools sort -@ ${task.cpus - 1} --write-index -o aligned.sorted.bam##idx##aligned.sorted.bam.bai
    samtools flagstats aligned.sorted.bam > ${meta.sample_name}_mapping_rates.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version 2>&1)"
        samtools: "\$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"
    END_VERSIONS
    """

    stub:
    """
    touch aligned.sorted.bam
    touch aligned.sorted.bam.bai
    touch ${meta.sample_name}_mapping_rates.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
        samtools: "stub"
    END_VERSIONS
    """
}
