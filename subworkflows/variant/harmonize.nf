// VARIANT_HARMONIZE — bring per-caller VCFs to a single shared schema so they
// can be compared and merged. BCFTOOLS_ANNOTATE_FILTER lifts the FILTER
// column into INFO/FILTER and BCFTOOLS_FILL_TAGS copies the site QUAL into
// INFO/VARQUAL; BCFTOOLS_INFO_TO_FORMAT then moves both INFOs into the
// per-caller FORMAT block so the sample-column-equals-caller scheme carries
// filter and quality per caller through the downstream merge. VAF is computed
// post-merge in BCFTOOLS_ANNOTATE_CONSENSUS, not here. Lofreq is site-only,
// so it still lacks a sample column at this point; we attach one by merging
// against a header-only TEMPLATE_VCF so downstream reheader+merge can proceed
// (the INFO→FORMAT lift is skipped for lofreq pending a separate fix).

include { BCFTOOLS_ANNOTATE_FILTER } from '../../modules/bcftools/annotate_filter/main'
include { BCFTOOLS_FILL_TAGS } from '../../modules/bcftools/fill_tags/main'
include { BCFTOOLS_INFO_TO_FORMAT } from '../../modules/bcftools/info_to_format/main'
include { TEMPLATE_VCF } from '../../modules/template_vcf/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_LOFREQ } from '../../modules/bcftools/merge/main'


workflow VARIANT_HARMONIZE {

    // input channels are tuples of (meta, caller, vcf.gz) — normalised per-caller VCFs.
    take:
        per_caller_vcf

    main:
        BCFTOOLS_ANNOTATE_FILTER(per_caller_vcf)
        BCFTOOLS_FILL_TAGS(BCFTOOLS_ANNOTATE_FILTER.out.vcf)

        BCFTOOLS_FILL_TAGS.out.vcf
            .branch { _meta, caller, _vcf ->
                lofreq: caller == 'lofreq'
                other:  true
            }
            .set { branched }

        BCFTOOLS_INFO_TO_FORMAT(branched.other)

        // Lofreq is site-only — attach a sample column via TEMPLATE_VCF + merge.
        TEMPLATE_VCF(branched.lofreq.map { meta, _caller, _vcf -> meta })

        branched.lofreq
            .combine(TEMPLATE_VCF.out.vcf, by: 0)
            .map { meta, _caller, lofreq_vcf, tpl_vcf, _tpl_tbi ->
                tuple(meta, [lofreq_vcf, tpl_vcf])
            }
            .set { lofreq_to_merge }

        BCFTOOLS_MERGE_LOFREQ(lofreq_to_merge)

        BCFTOOLS_MERGE_LOFREQ.out.vcf
            .map { meta, vcf -> tuple(meta, 'lofreq', vcf) }
            .set { lofreq_harmonized }

        harmonized = BCFTOOLS_INFO_TO_FORMAT.out.vcf.mix(lofreq_harmonized)

    // harmonized: (meta, caller, vcf.gz) — schema-aligned per-caller VCFs.
    emit:
        harmonized
}
