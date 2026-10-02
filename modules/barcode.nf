#!/usr/bin/env nextflow

process EXTRACT_BARCODE {
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    debug false

    input:
    tuple val(library_id), val(library_type), path(read2)
    val validation_done
    // this is just a signal to make sure the validation is done before this process starts
    val r2_lengths

    output:
    tuple path(read2), val(library_id), val(library_type), path("*_barcodes.parquet"), emit: barcode_extracted

    script:
    def sequence_length = r2_lengths[library_type.toLowerCase()]
    """
        if [ ${library_type} = DNA ]; then
            zcat ${read2} | awk '
            NR % 4 == 1 { 
                # Line 1 of each FASTQ record: @read_id
                split(substr(\$0, 2), name_parts, " ")
                read_id = name_parts[1];
            }
            NR % 4 == 2 {  # Line 2: sequence
            if (length(\$0) < '${sequence_length}') {
                    print "Error: Sequence length (" length(\$0) ") is less than required minimum ('${sequence_length}') for read " read_id > "/dev/stderr";
                    exit 1;
                }
                seq = substr(\$0, 1, ${sequence_length});
                c1 = substr(seq, 1, 8);
                c2 = substr(seq, 13, 8);
                c3 = substr(seq, 25, 8);
                c4 = substr(seq, 37, 8);
                printf "%s,%s,%s,%s,%s\\n", read_id, c1, c2, c3, c4;
            }
            '
        else
            zcat ${read2} | awk '
            NR % 4 == 1 { 
                # Line 1 of each FASTQ record: @read_id
                split(substr(\$0, 2), name_parts, " ")
                read_id = name_parts[1];
            }
            NR % 4 == 2 {  # Line 2: sequence
            if (length(\$0) < '${sequence_length}') {
                    print "Error: Sequence length (" length(\$0) ") is less than required minimum ('${sequence_length}') for read " read_id > "/dev/stderr";
                    exit 1;
                }
                seq = substr(\$0, 1, ${sequence_length});
                c1 = substr(seq, 1, 8);
                c2 = substr(seq, 13, 8);
                c3 = substr(seq, 25, 8);
                c4 = substr(seq, 37, 8);
                primer = substr(seq, 45, 18);
                umi = substr(seq, 63, 8);
                printf "%s,%s,%s,%s,%s,%s,%s\\n", read_id, c1, c2, c3, c4, primer, umi;
            }
            '
        fi > ${library_id}_${read2.baseName}.csv;
        convert_csv_to_parquet.py ${library_id}_${read2.baseName}.csv ${library_type} ${library_id};
        rm ${library_id}_${read2.baseName}.csv;
        """

    stub:
    """
    touch ${library_id}_${read2.baseName}_barcodes.parquet
    """
}

process CORRECT_BARCODE {
    label 'lazy_polars'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple path(read2), val(library_id), val(library_type), path(barcode_parquet_path), path(barcode_whitelist_path)

    output:

    tuple val(library_id), val(library_type), path("*_corrected_barcode_reads.parquet"), path("*_substitutions.parquet"), emit: corrected_barcodes
    tuple val(library_id), val(library_type), path("*_barcodes_mapping.parquet"), emit: barcode_mapping
    tuple val(library_id), path("metrics.parquet"), emit: metrics_parquet

    script:
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        correct_barcode.py ${barcode_parquet_path} ${barcode_whitelist_path} ${library_type} ${library_id}
        """

    stub:
    """
    mkdir -p ${library_id}
    touch ${library_id}_corrected_barcode_reads.parquet
    touch ${library_id}_barcodes_mapping.parquet
    touch metrics.parquet
    touch ${library_id}_substitutions.parquet
    """
}

process PREPARE_DEMULTIPLEXING {
    label 'small_job'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    path samplesheet_csv
    path barcode_whitelist_path
    val validation_done

    output:

    path ("sample_demultiplexing_sheet.csv"), emit: demultiplexing_csv
    path ("sample_sheet.csv"), emit: sample_sheet_csv

    script:
    """
        prepare_demultiplexing.py ${samplesheet_csv} ${barcode_whitelist_path} A
        """

    stub:
    """
    touch sample_demultiplexing_sheet.csv
    touch sample_sheet.csv
    """
}

process QC_BARCODE {
    label 'lazy_polars'
    tag "${library_id}"
    container "community.wave.seqera.io/library/uniflow-qc:c879a1b9c4c3978b"
    debug false

    input:
    tuple val(library_id), val(library_type), path(corrected_barcode_reads_parquet), path(mutations_parquet), path(barcode_whitelist_path)

    output:

    tuple val(library_id), path("*.png"), emit: barcode_qc

    script:
    """
        export NUMBA_CACHE_DIR="tmp/numba_cache"
        qc_barcode.py ${corrected_barcode_reads_parquet} ${mutations_parquet} ${barcode_whitelist_path} ${library_id} ${library_type}
        """

    stub:
    """
    touch ${library_id}_barcode_occurence_heatmap.png
    touch ${library_id}_barcode_upset_plot.png
    """
}
