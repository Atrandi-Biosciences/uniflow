include { DEMULTIPLEX_PARQUET ; DEMULTIPLEX_FASTQ ; CONCATENATE_FASTQ ; CONCATENATE_PARQUET } from '../modules/shared.nf'

include { PREPARE_DEMULTIPLEXING } from '../modules/barcode.nf'

workflow DEMULTIPLEXING {
    take:
    samples_ch
    corrected_barcodes_ch
    input_csv
    inputs_validated
    bc_long
    chemistry_json_ch

    main:


    PREPARE_DEMULTIPLEXING(input_csv, bc_long, inputs_validated)

    samples_ch
        | map { meta, read1, read2 ->
            [meta.library_id, meta.library_type, meta, read1, read2]
        }
        | combine(corrected_barcodes_ch, by: [0, 1])
        | combine(PREPARE_DEMULTIPLEXING.out.demultiplexing_csv)
        | DEMULTIPLEX_PARQUET


    DEMULTIPLEX_PARQUET.out.per_sample_reads
        | branch { meta, _read1, _read2, _per_sample_reads ->
            to_demultiplex: meta.demultiplexing_indices != null
            simple: true
        }
        | set { per_samples }

    per_samples.to_demultiplex
        | flatMap { meta, read1, read2, per_sample_reads ->
            [
                [meta, read1, per_sample_reads, "R1"],
                [meta, read2, per_sample_reads, "R2"],
            ]
        }
        | DEMULTIPLEX_FASTQ

    DEMULTIPLEX_FASTQ.out.demultiplexed_sample_fastq
        | set { demultiplexed_samples }

    per_samples.simple
        | map { meta, read1, read2, _per_sample_reads ->
            [meta, read1, read2]
        }
        | set { original_samples }

    original_samples.mix(
        demultiplexed_samples
            | groupTuple(by: 0)
            | map { meta, _strandedness, reads ->
                //use grep to identify R1 anf R2 files and assign read identities accordingly
                def read1 = reads.find { it.baseName.contains('_R1') || it.baseName.contains('read1') || it.baseName.contains('Read1') } ?: error("Could not identify read 1 file for sample ${meta.sample_name}")
                def read2 = reads.find { it.baseName.contains('_R2') || it.baseName.contains('read2') || it.baseName.contains('Read2') } ?: error("Could not identify read 2 file for sample ${meta.sample_name}")
                [meta, read1, read2]
            }
    )
        | set { pre_concatenation_samples }


    pre_concatenation_samples
        | map { meta, read1, read2 ->
            [meta - meta.subMap("library_id"), read1, read2]
        }
        | groupTuple(by: [0])
        | branch { meta, read1, read2 ->
            to_concatenate: read1.size() > 1 && read2.size() > 1
            no_concatenation: read1.size() == 1 && read2.size() == 1
        }
        | set { concatenation_branches }

    concatenation_branches.to_concatenate
        | CONCATENATE_FASTQ

    concatenation_branches.no_concatenation
        | map { meta, read1, read2 ->
            [meta, read1[0], read2[0]]
        }
        | mix(CONCATENATE_FASTQ.out.concatenated_fastq)
        | branch { meta, _read1, _read2 ->
            rna: meta.library_type == 'RNA'
            dna: meta.library_type == 'DNA'
        }
        | set { final_samples }


    DEMULTIPLEX_PARQUET.out.per_sample_barcode_mappings_parquet
        | map { meta, per_sample_mappings ->
            [meta - meta.subMap("library_id"), per_sample_mappings]
        }
        | groupTuple(by: [0])
        | branch { meta, per_sample_mappings ->
            to_concatenate: per_sample_mappings.size() > 1
            no_concatenation: per_sample_mappings.size() == 1
        }
        | set { concatenation_branches_parquet }

    concatenation_branches_parquet.to_concatenate
        | CONCATENATE_PARQUET

    concatenation_branches_parquet.no_concatenation
        | map { meta, per_sample_mappings ->
            [meta, per_sample_mappings[0]]
        }
        | mix(CONCATENATE_PARQUET.out.concatenated_parquet)
        | set { per_sample_barcode_mappings }

    emit:
    rna_libraries               = final_samples.rna
    dna_libraries               = final_samples.dna
    per_sample_barcode_mappings = per_sample_barcode_mappings
    metadata                    = PREPARE_DEMULTIPLEXING.out.sample_sheet_csv
}
