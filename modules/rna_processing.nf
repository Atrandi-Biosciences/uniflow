process ALIGN_GENE_EXPRESSION {

    //
    // This module executes STAR align quantification
    //
    label 'lazy_polars'
    container 'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7'
    tag "${meta.sample_name}"

    input:
    tuple val(meta), path(read1), path(read2)
    path rna_index
    path whitelist
    val protocol
    val star_feature
    val modality_name

    output:
    tuple val(meta), path('*d.out.bam'), emit: bam
    tuple val(meta), path('*.Solo.out'), emit: counts, optional: true
    tuple val(meta), path("*.Solo.out/Gene*/raw"), emit: raw_counts, optional: true
    tuple val(meta), path("*.Solo.out/Gene*/filtered"), emit: filtered_counts, optional: true
    tuple val(meta), path("*.Solo.out/Gene*/Summary.csv"), emit: solo_metrics_csv, optional: true
    tuple val(meta), path('*Log.final.out'), emit: log_final_txt
    tuple val(meta), path('*Log.out'), emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path "versions.yml", topic: versions
    tuple val(meta), path('*Aligned.sortedByCoord.out.bam'), optional: true, emit: bam_sorted
    tuple val(meta), path('*fastq.gz'), optional: true, emit: fastq
    tuple val(meta), path('*.tab'), optional: true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.sample_name}"

    // separate forward from reverse pairs
    def read_order = meta.library_type == 'RNA' ? "${read1} ${read2}" : "${read1} ${read2}"
    //def index_spec = meta.library_type == 'RNA' ? "${rna_index}" : "${dna_index}"

    """
    STAR \\
        --genomeDir ${rna_index} \\
        --readFilesIn ${read_order} \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        --soloCBwhitelist  ${whitelist} \\
        --soloType ${protocol} \\
        --soloCBposition 0_0_0_7 0_12_0_19 0_24_0_31 0_36_0_53 \\
        --soloUMIposition 0_62_0_69 \\
        --soloCBmatchWLtype EditDist_2 \\
        --soloFeatures ${star_feature} \\
        --outSAMtype BAM SortedByCoordinate \\
        --runDirPerm All_RWX \\
        --readFilesCommand zcat \\
        --outSAMattributes NH HI nM AS CR UR UB CB sS sQ sM GX GN \\
        --outSAMmultNmax 1 \\
        --limitBAMsortRAM ${task.memory.bytes} \\
        --outSAMunmapped Within

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    if [ -d ${prefix}.Solo.out ]; then
        # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "\$(STAR --version | sed -e "s/STAR_//g")"
    END_VERSIONS
    """

    stub:
    def prefix = "${modality_name}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "stub"
    END_VERSIONS
    """
}

process EXTRACT_GENE_EXPRESSION_BAM {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(bam)

    output:

    tuple val(meta), path("*_bam.parquet"), emit: bam_df_parquet
    tuple val(meta), path("*.bai"), emit: index_bai

    script:

    """
        extract_gene_expression_bam.py ${bam} ${meta.sample_name} ${meta.library_type} ${task.cpus}
        """

    stub:
    """
    touch ${meta.sample_name}_bam.parquet
    touch ${bam}.bai
    """
}

process GET_SATURATION_CURVES {
    label 'lazy_polars'
    tag "${meta.sample_name}"

    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(full_read_parquet), path(top_cells), val (internal_flag)

    output:

    tuple val(meta), path("*_saturation_metrics.csv"), emit: saturation_csv
    tuple val(meta.sample_name), path("*_saturation_mqc.png"), topic: per_sample_report
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    script:

    """
        get_saturation_curves.py ${full_read_parquet} ${top_cells} ${meta.sample_name} ${internal_flag}
        """

    stub:
    """
    touch ${meta.sample_name}_saturation_metrics.csv
    touch ${meta.sample_name}_saturation.png
    touch metrics.parquet
    """
}

process QC_AMPLICON_CARRYOVER {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug false

    input:
    tuple val(meta), path(dna_bam), path(rna_full_reads), path(top_cells_rna), path(reference_fasta)

    output:

    tuple val(meta), path("*.csv"), emit: metrics_csv
    tuple val(meta), path("*.png"), emit: carryover_amplicon_qc
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet

    script:

    """
        qc_amplicon_carryover.py ${dna_bam} ${rna_full_reads} ${top_cells_rna} ${reference_fasta} ${meta.sample_name}
        """

    stub:
    """
    touch ${meta.sample_name}_carryover.csv
    touch ${meta.sample_name}_carryover.png
    touch metrics.parquet
    """
}

process FILTER_GENE_EXPRESSION_CELLS {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug false

    input:

    tuple val(meta), path(full_read_parquet)

    output:
    tuple val(meta), path(full_read_parquet), path("*_top_cells.parquet"), emit: to_saturation, optional: true
    tuple val(meta), path("*_top_cells.parquet"), emit: top_cells_parquet, optional: true
    tuple val(meta), path("status.parquet"), optional: true, topic: "status"
    tuple val(meta), path("raw_counts.parquet"), emit: raw_counts_df_parquet, optional: true
    tuple val(meta), path("*.png"), emit: qc_plots, optional: true
    tuple val(meta.sample_name), path("*_knee_plot_mqc.png"), optional: true, topic: per_sample_report
    tuple val(meta), path("metrics.parquet"), emit: metrics_parquet
    tuple val(meta), path("raw_gene_counts.h5ad"), emit: raw_gene_counts_h5ad, optional: true

    script:
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        filter_gene_expression_cells.py ${full_read_parquet} ${meta.sample_name} ${meta.force_cells}
        """

    stub:
    """
    touch top_cells.parquet
    touch raw_counts.parquet
    touch stub_qc.png
    touch stub_metrics.parquet
    touch raw_gene_counts.h5ad
    touch status.parquet
    touch metrics.parquet
    """
}
