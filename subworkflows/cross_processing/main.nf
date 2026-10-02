include { COUNT_AMPLICON ; QC_AMPLICON ; EXTRACT_AMPLICON_BAM } from '../../modules/dna_processing.nf'
include { ALIGN_AMPLICON        } from '../../modules/minimap2/main'
include { MERGE_READS           } from '../../modules/shared.nf'
include { QC_AMPLICON_CARRYOVER } from '../../modules/rna_processing.nf'


workflow CROSS_PROCESSING {
    take:
    rna_libraries
    amplicon_fasta_file_ch
    rna_full_reads
    top_cells_parquet

    main:
    def modality_name = "cross_processing"
    ALIGN_AMPLICON(
        rna_libraries.map { meta, read1, _read2 -> [meta, read1] },
        amplicon_fasta_file_ch,
    )

    ALIGN_AMPLICON.out.bam
        | EXTRACT_AMPLICON_BAM

    EXTRACT_AMPLICON_BAM.out.dna_bam_parquet
        | map { meta, dna_bam ->
            [meta - meta.subMap("library_type", "library_id"), dna_bam]
        }
        | combine(
            rna_full_reads | map { meta, full_reads ->
                [meta - meta.subMap("library_type", "library_id"), full_reads]
            },
            by: 0
        )
        | combine(
            top_cells_parquet | map { meta, top_cells ->
                [meta - meta.subMap("library_type", "library_id"), top_cells]
            },
            by: 0
        )
        | combine(amplicon_fasta_file_ch)
        | QC_AMPLICON_CARRYOVER

    emit:
    metrics         = channel.empty().groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    internal_metrics = channel.empty().mix(ALIGN_AMPLICON.out.mapping_rates).mix(QC_AMPLICON_CARRYOVER.out.metrics_csv).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    qc              = channel.empty().mix(QC_AMPLICON_CARRYOVER.out.carryover_amplicon_qc).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
    metrics_parquet = QC_AMPLICON_CARRYOVER.out.metrics_parquet
    dev             = channel.empty()
    out             = channel.empty()
}
