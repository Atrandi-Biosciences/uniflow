// VARIANT_CALLING — bulk variant-calling subworkflow.

include { SAMTOOLS_FAIDX                    } from '../../modules/samtools/faidx/main'
include { BCFTOOLS_FILTER                   } from '../../modules/bcftools/filter/main'
include { BCFTOOLS_NORM                     } from '../../modules/bcftools/norm/main'
include { BCFTOOLS_STATS                    } from '../../modules/bcftools/stats/main'
include { BCFTOOLS_CALL                     } from '../../modules/bcftools/call/main'
include { CONCAT_VCF                        } from '../../modules/bcftools/concat/main'
include { FREEBAYES                         } from '../../modules/freebayes/main'
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../../modules/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER             } from '../../modules/gatk4/haplotypecaller/main'
include { GATK4_UPDATEVCFSEQUENCEDICTIONARY } from '../../modules/gatk4/updatevcfsequencedictionary/main'
include { LOFREQ_CALL                       } from '../../modules/lofreq/call/main'
include { BEDTOOLS_MAKEWINDOWS              } from '../../modules/bedtools/makewindows/main'
include { TABIX                             } from '../../modules/tabix/tabix/main'
include { VARDICT                           } from '../../modules/vardict/main'


workflow VARIANT_CALLING {
    take:
    bam_with_index
    amplicon_fasta

    main:
    def callers = params.variant_callers.tokenize(',')
    def unsupported = callers - ['freebayes', 'lofreq', 'bcftools', 'gatkhc', 'vardict']
    if (unsupported) {
        error("variant_callers contains ${unsupported}")
    }

    SAMTOOLS_FAIDX(amplicon_fasta)

    // Amplicon BED: built unconditionally now, since FreeBayes + GATK scatter over it.
    // VarDict needs the whole BED, and it is published next to the BAM.
    // TODO: if amplicon_bed is provided, we should validate that it matches the fasta
    def amplicon_bed_ch = null
    if (params.amplicon_bed) {
        amplicon_bed_ch = channel.value(file(params.amplicon_bed, checkIfExists: true))
    }
    else {
        BEDTOOLS_MAKEWINDOWS(SAMTOOLS_FAIDX.out.fai)
        amplicon_bed_ch = BEDTOOLS_MAKEWINDOWS.out.bed
    }

    // 1. Scatter: per-amplicon intervals for the scatter-gather callers (FreeBayes, GATK HC):
    // one interval per BED line = (interval_id, bed_line_text)
    def scatter_callers = ['freebayes', 'gatkhc']
    def scatter_input = null
    if (scatter_callers.any { it in callers }) {
        def intervals = amplicon_bed_ch
            .splitText()
            .map { it.trim() }
            .filter { it }
            .map { line ->
                def cols = line.split('\t')
                // interval_id we create here must be unique per BED line (not just per contig)
                // so per-interval VCF filenames never collide
                tuple("${cols[0]}_${cols[1]}_${cols[2]}", line)
            }
        scatter_input = bam_with_index.combine(intervals)
    }
    def caller_vcfs = channel.empty()
    def needs_dict_update = channel.empty()
    // scattered_vcfs stays a build-time list: its emptiness gates whether
    // CONCAT_VCF is instantiated at all (a channel would always be truthy).
    def scattered_vcfs = []

    if ('freebayes' in callers) {
        FREEBAYES(
            scatter_input,
            SAMTOOLS_FAIDX.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
        )
        scattered_vcfs << FREEBAYES.out.vcf
    }
    if ('lofreq' in callers) {
        LOFREQ_CALL(
            bam_with_index,
            SAMTOOLS_FAIDX.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
        )
        needs_dict_update = needs_dict_update.mix(LOFREQ_CALL.out.vcf)
    }
    if ('bcftools' in callers) {
        BCFTOOLS_CALL(
            bam_with_index,
            SAMTOOLS_FAIDX.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
        )
        caller_vcfs = caller_vcfs.mix(BCFTOOLS_CALL.out.vcf)
    }
    if ('gatkhc' in callers) {
        GATK4_CREATESEQUENCEDICTIONARY(SAMTOOLS_FAIDX.out.fasta)
        GATK4_HAPLOTYPECALLER(
            scatter_input,
            SAMTOOLS_FAIDX.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
            GATK4_CREATESEQUENCEDICTIONARY.out.dict,
        )
        scattered_vcfs << GATK4_HAPLOTYPECALLER.out.vcf
    }
    // 2. Gather: group the per-interval VCFs by (sample, caller)
    // A sort is added to the grouped list to turn it into a canonical (by-filename) order.
    // This is purely cosmetic, so the CONCAT_VCF task hash is stable and can produce the same
    // result with -resume. Output hash will always be ordered!
    if (scattered_vcfs) {
        Channel.empty().mix(*scattered_vcfs)
            | groupTuple(by: [0, 1])
            | map { meta, caller, vcfs -> tuple(meta, caller, vcfs.toSorted { it.name }) }
            | set { scattered_grouped }
        CONCAT_VCF(scattered_grouped)
        caller_vcfs = caller_vcfs.mix(CONCAT_VCF.out.vcf)
    }
    if ('vardict' in callers) {
        VARDICT(
            bam_with_index,
            SAMTOOLS_FAIDX.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
            amplicon_bed_ch,
        )
        needs_dict_update = needs_dict_update.mix(VARDICT.out.vcf)
    }

    GATK4_UPDATEVCFSEQUENCEDICTIONARY(needs_dict_update.combine(bam_with_index, by: 0))
    caller_vcfs = caller_vcfs.mix(GATK4_UPDATEVCFSEQUENCEDICTIONARY.out.vcf)

    BCFTOOLS_NORM(
        caller_vcfs,
        SAMTOOLS_FAIDX.out.fasta,
        SAMTOOLS_FAIDX.out.fai,
    )

    BCFTOOLS_FILTER(BCFTOOLS_NORM.out.vcf)

    TABIX(BCFTOOLS_FILTER.out.vcf)
    BCFTOOLS_STATS(BCFTOOLS_FILTER.out.vcf)

    emit:
    vcf       = BCFTOOLS_FILTER.out.vcf
    tbi       = TABIX.out.tbi
    stats     = BCFTOOLS_STATS.out.stats
    fasta     = SAMTOOLS_FAIDX.out.fasta
    fai       = SAMTOOLS_FAIDX.out.fai
    out_stats = BCFTOOLS_STATS.out.out_stats
}
