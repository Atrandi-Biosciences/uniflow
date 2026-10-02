process DEMULTIPLEX_PARQUET {
    label 'lazy_polars'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(library_id), val(library_type), val(meta), path(read1), path(read2), path(barcode_mappings_parquet), path(demultiplexing_sheet)

    output:
    tuple val(meta), path("*_barcode_mappings.parquet"), emit: per_sample_barcode_mappings_parquet
    tuple val(meta), path(read1), path(read2), path("*_reads.txt"), emit: per_sample_reads

    script:
    """
        demultiplex_parquet.py ${library_id} ${library_type} ${barcode_mappings_parquet} ${demultiplexing_sheet} ${read1} ${read2} ${meta.sample_name}
        """

    stub:
    """
    touch ${meta.sample_name}_corrected_barcode_reads.parquet
    touch ${meta.sample_name}_barcode_mappings.parquet
    touch ${meta.sample_name}_reads.txt
    """
}

process CONCATENATE_PARQUET {
    label 'cpu_low_mem'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(parquet_files, stageAs: "?/*")

    output:
    tuple val(meta), path("*_barcode_mappings.parquet"), emit: concatenated_parquet

    script:
    // Build a safe, space-separated list of quoted paths
    """
    concatenate_parquet.py ${meta.sample_name} ${meta.library_type} "${parquet_files}"
    """
}

process DEMULTIPLEX_FASTQ {
    label 'demultiplexing_fastq'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/seqtk_pigz:aa99a20f06d8e9a8"
    debug false

    input:
    tuple val(meta), path(fastq), path(per_sample_reads), val(read_identity)

    output:
    tuple val(meta), val(read_identity), path("*.fastq.gz"), emit: demultiplexed_sample_fastq

    script:
    def output_name = "${meta.sample_name}_${meta.library_type}_${read_identity}.fastq.gz"
    """
        set -o pipefail
        seqtk subseq ${fastq} ${per_sample_reads} | pigz -p ${task.cpus} > ${output_name}
        """

    stub:
    def suffix = fastq.baseName.contains('_R1') ? 'R1' : fastq.baseName.contains('_R2') ? 'R2' : error('Fastq filename must contain "R1" or "R2" in the basename')
    def output_name = "${meta.sample_name}_${meta.library_type}_${suffix}.fastq.gz"
    """
    touch ${output_name}
    """
}

process CONCATENATE_FASTQ {
    label 'small_job'
    tag "${meta.sample_name}"
    debug false

    input:
    tuple val(meta), path(read1_files, stageAs: "?/*"), path(read2_files, stageAs: "?/*")

    output:
    tuple val(meta), path("*_R1.fastq.gz"), path("*_R2.fastq.gz"), emit: concatenated_fastq

    script:

    """
        cat ${read1_files.join(' ')} > ${meta.sample_name}_${meta.library_type}_R1.fastq.gz
        cat ${read2_files.join(' ')} > ${meta.sample_name}_${meta.library_type}_R2.fastq.gz
    """
}

process MERGE_H5AD {
    label 'small_job'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    debug false

    input:
    // Accept a variadic list of .h5ad paths as a single list-valued path input
    tuple val(meta), path(sample_sheet), val(experiment_id), path(h5ad_paths)

    output:
    tuple val(meta), path("counts.h5mu"), emit: counts_h5ad
    tuple val(meta), path("counts.zarr"), emit: counts_zarr

    script:
    // Build a safe, space-separated list of quoted paths
    """
    
    merge_h5ad.py ${meta.sample_name} ${sample_sheet} ${experiment_id} ${h5ad_paths}
    """

    stub:
    """
    touch counts.h5mu
    touch counts.zarr
    """
}


process MERGE_METRICS {
    label 'small_job'
    tag "${source_id}"
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    debug false

    input:
    tuple val(source_id), path(metrics_paths, stageAs: "?/*")

    output:
    tuple val(source_id), path("metrics.csv"), emit: metrics_csv

    script:
    """
        merge_metrics.py ${source_id} ${metrics_paths}
    """

    stub:
    """
    touch metrics.csv
    """
}

process MERGE_STATUS {
    label 'small_job'
    tag "${source_id}"
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    debug false

    input:
    tuple val(source_id), path(status_parquet_paths, stageAs: "?/*")

    output:
    tuple val(source_id), path("status.csv"), emit: status_parquet

    script:
    """
        merge_status.py ${source_id} ${status_parquet_paths}
    """

    stub:
    """
    touch status.csv
    """
}

// Take the barcode parquet and the bam extracted parquet to merge the reads and provide a full raw dataframe of R1 and R2 reads mapped by read_name
process MERGE_READS {
    label 'lazy_polars'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(bam_parquet), path(barcodes_mapping_parquet)

    output:
    tuple val(meta), path("full_read_*.parquet"), emit: merged_reads_parquet

    script:
    """
        merge_reads.py ${bam_parquet} ${barcodes_mapping_parquet} ${meta.sample_name} ${meta.library_type}
        """

    stub:
    """
    touch full_read_${meta.library_type}.parquet
    """
}

process TRIM_READS {
    label 'cpu_low_mem'
    tag "${meta.sample_name}"
    container "quay.io/biocontainers/cutadapt:1.18--py35_0"

    debug false

    input:
    tuple val(meta), path(read1), path(read2)
    val trimming_length

    output:
    tuple val(meta), path("*R1.trimmed.fastq.gz"), path("*R2.trimmed.fastq.gz"), emit: trimmed_fastq

    script:

    //def TSO_SEQ = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    """
        cutadapt  --length ${trimming_length} --cores ${task.cpus} --output ${meta.sample_name}_R1.trimmed.fastq.gz --paired-output ${meta.sample_name}_R2.trimmed.fastq.gz ${read1} ${read2} > ${meta.sample_name}_cutadapt_tso.log
        """

    stub:
    """
    touch ${meta.sample_name}_R1.trimmed.fastq.gz
    touch ${meta.sample_name}_R2.trimmed.fastq.gz
    """
}
