include { VALIDATE_AMPLICON_FASTA } from '../../modules/validation.nf'


workflow FASTA_VALIDATION {
    take:
    amplicon_fasta_file_ch

    main:
    VALIDATE_AMPLICON_FASTA(amplicon_fasta_file_ch)

    emit:
    done       = VALIDATE_AMPLICON_FASTA.out.report.map { true }.ifEmpty(true)
    similarity = VALIDATE_AMPLICON_FASTA.out.similarity
    report     = VALIDATE_AMPLICON_FASTA.out.report
}
