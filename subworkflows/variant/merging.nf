// VARIANT_MERGING — collapse per-caller harmonised VCFs into a single
// per-sample VCF and annotate it with the set of callers that detected
// each site (INFO/CALLERS, INFO/NCALLERS).

include { BCFTOOLS_REHEADER_SAMPLE } from '../../modules/bcftools/reheader_sample/main'
include { BCFTOOLS_MERGE } from '../../modules/bcftools/merge/main'
include { BCFTOOLS_ANNOTATE_CONSENSUS } from '../../modules/bcftools/annotate_consensus/main'
include { TABIX } from '../../modules/tabix/tabix/main'


workflow VARIANT_MERGING {

    // input channels are tuples of (meta, caller, vcf.gz) — harmonised per-caller VCFs.
    take:
        harmonized

    main:
        // Rename each sample column to the caller string, group by sample,
        // then bcftools merge into one VCF per sample.
        BCFTOOLS_REHEADER_SAMPLE(harmonized)

        BCFTOOLS_REHEADER_SAMPLE.out.vcf
            .map { meta, _caller, vcf -> tuple(meta, vcf) }
            .groupTuple(by: 0)
            .set { per_sample_reheaded }

        BCFTOOLS_MERGE(per_sample_reheaded)

        BCFTOOLS_ANNOTATE_CONSENSUS(BCFTOOLS_MERGE.out.vcf)
        
        BCFTOOLS_MERGE.out.vcf
            .map { meta, vcf -> tuple(meta, 'merged', vcf) }
            .mix(BCFTOOLS_ANNOTATE_CONSENSUS.out.vcf
                .map { meta, vcf -> tuple(meta, 'consensus', vcf) })
            .set { to_index }

        TABIX(to_index)

        TABIX.out.tbi
            .branch { _meta, label, _tbi ->
                merged:    label == 'merged'
                consensus: label == 'consensus'
            }
            .set { tbi_branched }

    // merged_vcf:    (meta, vcf.gz)       per-sample bcftools merge output
    // merged_tbi:    (meta, vcf.gz.tbi)
    // consensus_vcf: (meta, vcf.gz)       merge_vcf with INFO/CALLERS, INFO/NCALLERS
    // consensus_tbi: (meta, vcf.gz.tbi)
    emit:
        merged_vcf = BCFTOOLS_MERGE.out.vcf
        merged_tbi = tbi_branched.merged.map { meta, _label, tbi -> tuple(meta, tbi) }
        consensus_vcf = BCFTOOLS_ANNOTATE_CONSENSUS.out.vcf
        consensus_tbi = tbi_branched.consensus.map { meta, _label, tbi -> tuple(meta, tbi) }
}
