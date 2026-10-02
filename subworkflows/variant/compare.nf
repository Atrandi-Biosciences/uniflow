// VARIANT_COMPARE - cross caller concordance QC. Takes the harmonised
// per-caller VCFs for each sample and emits a comparison report.

include { CALLER_COMPARE_QC } from '../../modules/qc/caller_compare/main'


workflow VARIANT_COMPARE {

    // harmonized:    (meta, caller, vcf.gz) — harmonised per-caller VCFs
    // amplicon_fasta: path
    // amplicon_fai:   path
    take:
        harmonized
        amplicon_fasta
        amplicon_fai

    main:
        harmonized
            .map { meta, caller, vcf -> tuple(meta, caller, vcf) }
            .groupTuple(by: 0)
            .set { per_sample_harmonized }

        CALLER_COMPARE_QC(
            per_sample_harmonized,
            amplicon_fasta,
            amplicon_fai,
        )

    // qc_report: (meta, dir)
    emit:
        qc_report = CALLER_COMPARE_QC.out.report
}
