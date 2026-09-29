from dataclasses import dataclass
from typing import Final

import numpy as np
import polars as pl
import pysam
from pysam.libcalignmentfile import IteratorColumnRegion
from scipy.special import betaln, gammaln, xlogy

from lib.common_const import BARCODE, READ_NAME, READS

# Per-cell attribution parameters (defaults only; the live values are Nextflow
# params passed as named flags to count_variant*.py. see nextflow.config /
# conf/variant.config). No hard drops remain in the per-cell layer except the
# pooled site gate; the former coverage / VAF checks are now soft FILTER
# labels (lowDP / lowVAF), nothing is deleted.
#
#   DEFAULT_MIN_SITE_POOLED_COVERAGE : minimum reads spanning a candidate site,
#       summed over all cells, for the per-cell pileup to attribute it (was the
#       hardcoded MIN_COVERAGE). The only remaining hard drop; a discovery-scope
#       guard, not a per-cell filter.
#
#   DEFAULT_ERROR_EPS / DEFAULT_ERROR_RHO / DEFAULT_MIN_GQ : the genotype
#       likelihood parameters. See "Note on the genotype" below; the defaults
#       are the values measured on the 10M run, and every run is expected to
#       measure its own.
#
#   DEFAULT_FLAG_MIN_DP : the depth below which a per-allele row is labelled
#       lowDP. It is also the cell-support callable-coverage denominator. It no
#       longer decides whether a genotype is callable; GQ does.
#
#   DEFAULT_FLAG_MIN_VAF : the per-allele lowVAF label in raw_variants and the
#       per-cell "vote" threshold for cell support. Deliberately NOT part of the
#       genotype FILTER: poor ALT evidence is a property of the allele, not a
#       reason to distrust the cell's call.
#
#   DEFAULT_FLAG_MIN_BQ / DEFAULT_FLAG_MIN_MQ : the per-allele lowBQ / lowMQ
#       labels in raw_variants, applied to the reported mean_bq / mean_mq of the
#       reads supporting that allele. Like lowVAF these are observation-quality
#       annotations and are deliberately NOT part of the genotype FILTER (see
#       genotype_filter_expr). A null stat is "unknown", never "low": indels
#       have no anchor base quality, and an all-reference cell-site has no
#       chosen ALT to take statistics from.
#       NOTE: 20.0 was calibrated against the OLD mean-of-Phred aggregate, which
#       Jensen inflated by a median of ~11 Phred on the dominant allele and ~2.5
#       on minor alleles (measured, 10M run). mean_bq_expr now reports the honest
#       value, so the same 20.0 is a materially stricter cut than it was and is
#       due a recalibration alongside the GT/GQ work.
#
# Note on the genotype:
# The genotype is the argmin of three Phred-scaled beta-binomial likelihoods, not a
# VAF band. The inherited 0.2 / 0.8 rule is gone: inverting the likelihood exactly,
# a hom-ref/het boundary at 0.2 encodes an error rate of ~3.6%, roughly 80x above
# the ~4.5e-4 this assay actually shows, and nobody chose that rate deliberately.
# GQ replaces depth as the callability gate,
# which is what Mission Bio Tapestri does (min_gq = 30) and what we could not do
# before, since we produced no per-cell genotype quality.

DEFAULT_MIN_SITE_POOLED_COVERAGE: Final[int] = 11
DEFAULT_FLAG_MIN_DP: Final[int] = 5
DEFAULT_FLAG_MIN_VAF: Final[float] = 0.01
DEFAULT_FLAG_MIN_BQ: Final[float] = 20.0
DEFAULT_FLAG_MIN_MQ: Final[float] = 30.0

# Genotype-likelihood parameters:
#
# eps and rho are MEASUREMENTS, not constants. These two defaults are the values
# measured on the 10M-read run and are seeds for a run that has not measured its
# own yet, thus we expects each run to re-estimate them. eps is a floor: it is
# the panel-wide machinery error, and the per-substitution-type table spans 23x
# around it.
#
# theta_het is fixed at 0.5 and is deliberately not a parameter. Measured 0.4977
# pooled and 0.469-0.548 per amplicon, and the measured cost of
# assuming perfect balance is 0.7 phred of median GQ. Revisit only if an amplicon
# drifts below ~0.35, which nothing on this run does.
#
# DEFAULT_MIN_GQ is a placeholder, not a calibrated value. percell_min_gq is a
# function of the run's cell count via a false-carrier limit:
# GQ >= 30 costs ~1-3 false carriers panel-wide at 20k cells but 10-40 at 300k.
# 30 matches Tapestri's min_gq and is where Stage D starts, not where it ends.
DEFAULT_ERROR_EPS: Final[float] = 4.5e-4
DEFAULT_ERROR_RHO: Final[float] = 0.091
DEFAULT_MIN_GQ: Final[float] = 30.0
THETA_HET: Final[float] = 0.5

# eps is per observation, not just per panel: the hom-ref arm's expected
# ALT fraction is raised to the ALT-supporting reads' own measured error
# probability whenever base quality says they are junk.
# A wrong base call is one of the three non-reference bases, so a read whose
# error probability is p supports any one specific ALT with probability p / 3.
# That factor is what converts the instrument's "this base is wrong" into the
# "this read shows this allele" the likelihood's theta actually means.
ALT_BASES_PER_SITE: Final[float] = 3.0

# VCF convention: GQ is capped at 99, so a saturated call reports 99 rather than
# the several-hundred-phred number the likelihood ratio actually reaches.
GQ_CAP: Final[float] = 99.0


CHROM: Final[str] = "chrom"
REF: Final[str] = "REF"
ALT: Final[str] = "ALT"
POS: Final[str] = "POS"
QUAL: Final[str] = "QUAL"
VAR_FEATURE: Final[str] = "var_feature"
FEATURE: Final[str] = "feature"
SNP_FEATURE: Final[str] = "snp_feature"
BASES: Final[str] = "bases"
DP: Final[str] = "DP"
READ_BASE: Final[str] = "read_base"
REF_BASE: Final[str] = "ref_base"
VARIANT_NAME: Final[str] = "variant_name"

# Dense cell key used inside the pileup loop. See build_read_cell_lut
CELL_INDEX: Final[str] = "cell_index"

# Downstream-friendly raw_variants.parquet columns
TOTAL_READS: Final[str] = "total_reads"
POS_LOCAL: Final[str] = "pos_local_0based"
REF_OUT: Final[str] = "ref"
ALT_OUT: Final[str] = "alt"
REF_READS: Final[str] = "ref_reads"
ALT_READS: Final[str] = "alt_reads"
VAF: Final[str] = "vaf"
IS_NO_CALL: Final[str] = "is_no_call"
N_CALLERS: Final[str] = "n_callers"
CALLERS: Final[str] = "callers"
ZYGOSITY: Final[str] = "zygosity"

# Genotype-likelihood output columns. PL is VCF-standard:
# -10*log10 P(k|n,G) per diploid state, shifted so the smallest is 0, so PL is a
# log-odds against the called genotype and the called state always reads 0. They
# stay Float64 rather than the VCF integer, since parquet costs nothing to keep
# the precision and a consumer that wants integers can round.
GT: Final[str] = "gt"
GQ: Final[str] = "gq"
PL_HOM_REF: Final[str] = "pl_homref"
PL_HET: Final[str] = "pl_het"
PL_HOM_ALT: Final[str] = "pl_homalt"
OUTLIER_LOD: Final[str] = "outlier_lod"

# VCF genotype strings, in the PL order the VCF spec fixes (0/0, 0/1, 1/1).
GT_HOM_REF: Final[str] = "0/0"
GT_HET: Final[str] = "0/1"
GT_HOM_ALT: Final[str] = "1/1"
GT_LABELS: Final[tuple[str, str, str]] = (GT_HOM_REF, GT_HET, GT_HOM_ALT)

# Per-cell quality statistics carried out of the pileup and reported (never used
# as a gate here; base-quality filtering belongs to the callers upstream). BASE_QUAL
# / MAP_QUAL are per-read; MEAN_BQ / MEAN_MQ are the per-cell aggregates. FILTER is
# the soft label: PASS / lowDP / lowVAF on the per-allele rows, PASS / lowGQ on
# the genotype rows; nothing is dropped for it either way.
BASE_QUAL: Final[str] = "base_qual"
MAP_QUAL: Final[str] = "map_qual"
MEAN_BQ: Final[str] = "mean_bq"
MEAN_MQ: Final[str] = "mean_mq"
FILTER: Final[str] = "filter"

# Soft-FILTER label values.
FILTER_PASS: Final[str] = "PASS"
FILTER_LOW_DP: Final[str] = "lowDP"
FILTER_LOW_VAF: Final[str] = "lowVAF"
FILTER_LOW_BQ: Final[str] = "lowBQ"
FILTER_LOW_MQ: Final[str] = "lowMQ"
FILTER_LOW_GQ: Final[str] = "lowGQ"
LOWQUAL_NAME: Final[str] = "lowqual"

# Genotype label values, the historic strings kept as the h5ad obs vocabulary.
# homref is a called state, not a missing one: a cell whose reads make 0/0 the
# most likely genotype at a usable GQ is confidently reference, which is a
# different statement from lowqual (not enough evidence to call anything).
ZYG_HOM_ALT: Final[str] = "hom"
ZYG_HOM_REF: Final[str] = "homref"
ZYG_HET: Final[str] = "het"
HET_NAME: Final[str] = "hete"


def mean_bq_expr(bq_col: str):
    """Aggregate per-read Phred base qualities into one Phred number, correctly.

    A Phred score is -10*log10(p_error), so averaging Phred scores averages
    logarithms: by Jensen's inequality that always OVERSTATES the quality of a
    mixed group, and it overstates it badly whenever the group is heterogeneous.
    Measured on this run's own BAM, mean(Phred) sits a median of 11 Phred above
    the honest value for a cell-sized group of reads on the dominant allele, and
    ~2.5 Phred above it for the thin minor-allele groups.

    The right aggregate is the mean ERROR PROBABILITY, converted back to Phred so
    the number keeps its units and stays comparable with flag_min_bq:

        mean_bq = -10 * log10( mean( 10^(-BQ/10) ) )

    Nulls are skipped by mean(), so a group whose reads all carry QUAL "*" still
    aggregates to null and stays "unknown" rather than "low" (see
    soft_filter_expr).
    """
    return (-10.0 * (pl.lit(10.0).pow(-pl.col(bq_col) / 10.0)).mean().log10()).alias(
        MEAN_BQ
    )


def soft_filter_expr(
    dp_col: str,
    vaf_col: str,
    flag_min_dp: int,
    flag_min_vaf: float,
    bq_col: str | None = None,
    mq_col: str | None = None,
    flag_min_bq: float = DEFAULT_FLAG_MIN_BQ,
    flag_min_mq: float = DEFAULT_FLAG_MIN_MQ,
):
    """The soft per-allele FILTER label (VCF-style, semicolon-joined reasons).

    PASS when the row clears every threshold, otherwise the semicolon-joined
    list of the reasons that fired, in a fixed order: lowDP, lowVAF, lowBQ,
    lowMQ (e.g. lowDP;lowMQ). Never used to drop a row -- it only annotates it.

    Depth and VAF treat a null as below the threshold: an all-reference cell has
    no ALT evidence, and saying so is the point. The quality statistics do the
    opposite and treat a null as *unknown*, so it never fires: indels carry no
    anchor base quality (mean_bq is null by construction) and a row with no
    chosen ALT has no reads to take statistics from. Labelling those lowBQ
    would report a measurement that was never made.

    Pass bq_col / mq_col to enable the quality reasons; omit either one when
    the frame does not carry it. This is the label for an *allele observation*
    (raw_variants.parquet). The genotype tables use genotype_filter_expr
    instead, since thin or noisy ALT evidence is not a defect of the call.
    """
    reasons: list[tuple[pl.Expr, str]] = [
        (pl.col(dp_col).fill_null(0) < flag_min_dp, FILTER_LOW_DP),
        (pl.col(vaf_col).fill_null(0.0) < flag_min_vaf, FILTER_LOW_VAF),
    ]
    if bq_col is not None:
        reasons.append((pl.col(bq_col) < flag_min_bq, FILTER_LOW_BQ))
    if mq_col is not None:
        reasons.append((pl.col(mq_col) < flag_min_mq, FILTER_LOW_MQ))

    parts = [
        pl.when(cond).then(pl.lit(label)).otherwise(pl.lit(None, dtype=pl.String))
        for cond, label in reasons
    ]
    joined = pl.concat_str(parts, separator=";", ignore_nulls=True)
    return (
        pl.when(joined.str.len_chars() == 0).then(pl.lit(FILTER_PASS)).otherwise(joined)
    )


def _beta_binomial_log_pmf(k, n, theta, rho: float):
    """Natural-log P(k ALT of n | expected ALT fraction theta) without C(n,k).

    Beta-binomial with alpha = theta*s, beta = (1-theta)*s, s = (1-rho)/rho, so
    rho is the intra-class correlation and Var(k) = n*theta*(1-theta)*[1+(n-1)*rho].

    C(n,k) is dropped because it is identical for all three diploid states at the
    same (k, n) and cancels out of every PL difference. It does not cancel against
    the uniform-theta outlier, which is why outlier_lod adds it back.

    rho = 0 is the binomial as an exact special case (s would be infinite, so).
    xlogy keeps 0*log(0) at 0, which matters at eps = 0 with no ALT reads.
    """
    if rho <= 0.0:
        return xlogy(k, theta) + xlogy(n - k, 1.0 - theta)
    s = (1.0 - rho) / rho
    alpha = theta * s
    beta = (1.0 - theta) * s
    return betaln(k + alpha, n - k + beta) - betaln(alpha, beta)


def genotype_likelihood_arrays(
    ref_reads: np.ndarray,
    alt_reads: np.ndarray,
    eps: float,
    rho: float,
    theta_het: float = THETA_HET,
    alt_base_error: np.ndarray | None = None,
):
    """PL / GQ / GT / outlier_lod for arrays of (ref_reads, alt_reads).

    n = ref_reads + alt_reads, not total_reads: a third allele is evidence about
    neither the REF nor the called ALT, and putting it in the denominator biases
    every call toward hom-ref; and we want to avoid that.

    alt_base_error is the per-observation error probablity of the ALT-supporting
    reads (10^(-mean_bq/10), not a Phred score), and it makes the hom-ref arm
    quality-awar. Without it a single junk base is enough to sink a
    well-covered reference cell. This is because our flat error eps set at 4.5e4.
    which means we trust a Q9 base ~280x more than the instrument does.

    Pass None (indels, which carry no anchor base quality) to get the flat-eps
    behaviour exactly.

    Some more notes:

    - This will only ever raise the hom-ref expectation. the eps stays a threshold because it is
    the panel-wide assay error rate and it is for sequencing, PCR, ambient and alignment error
    that base quality alone cannot see. And because a Q40 read's 1e-4 is below what the
    assay demonstrably achieves. That one-sidedness is what keeps the false carrier-cell limit
    valid as an upper bound (the rate is monotone decreasing in the calling eps), so the report
    needs no compensating change.

    - The hom-alt arm is deliberately left on the flat eps. Its error term is the
    REF-supporting reads' quality, which the genotype row does not carry. Using
    the ALT reads' quality there would apply a statistic to the wrong arm.

    - The three diploid states differ only in their expected ALT fraction. The eps
    (raised per observation) for hom-ref, theta_het for het, 1 - eps for hom-alt.
    Returns, in order, (pl_homref, pl_het, pl_homalt, gq, gt_index, outlier_lod),
    where the PLs are shifted so the smallest is 0, gq is the second-smallest PL
    capped at GQ_CAP, and gt_index indexes GT_LABELS.

    Ties go to the lower index (hom-ref, then het), which is what makes an
    evidence-free row (n = 0) come out 0/0 at GQ 0 rather than at random. This means
    all three likelihoods are exactly 1 there, and GQ 0 is the correct report.

    outlier_lod is log10 P(k | n, theta ~ Uniform) - log10 P(k | n, best diploid),
    essentially a way to say "this cell is not any diploid genotype" (doublet, ambient,
    early PCR over amplification, subclonal-within-cell, etc). Marginalising theta over
    Uniform(0,1) gives exactly 1/(n+1) independent of k, which is why it needs no
    parameter here. But it is a full probability, so the diploid side must carry its
    C(n,k) back. It is basically a score and it touches neither GT nor GQ.
    """
    k = np.asarray(alt_reads, dtype=np.float64)
    n = k + np.asarray(ref_reads, dtype=np.float64)

    theta_homref = eps
    if alt_base_error is not None:
        # A null mean_bq (no ALT allele in this cell, or an indel) comes out as NaN
        # and means "not measured", so it must fall back to the floor rather than
        # propagate. So here nan_to_num is used to assign it to a 0 error probability.
        per_obs = np.nan_to_num(np.asarray(alt_base_error, dtype=np.float64), nan=0.0)
        theta_homref = np.maximum(eps, per_obs / ALT_BASES_PER_SITE)

    log_p = np.vstack(
        [
            _beta_binomial_log_pmf(k, n, theta_homref, rho),
            _beta_binomial_log_pmf(k, n, theta_het, rho),
            _beta_binomial_log_pmf(k, n, 1.0 - eps, rho),
        ]
    )

    ln10 = np.log(10.0)
    phred = (-10.0 / ln10) * log_p
    phred -= phred.min(axis=0)
    gt_index = np.argmin(phred, axis=0)
    gq = np.minimum(np.sort(phred, axis=0)[1], GQ_CAP)

    log10_choose = (gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)) / ln10
    log10_best_diploid = log10_choose + log_p.max(axis=0) / ln10
    outlier_lod = -np.log10(n + 1.0) - log10_best_diploid

    return phred[0], phred[1], phred[2], gq, gt_index, outlier_lod


GENOTYPE_LIKELIHOOD_SCHEMA: Final[dict] = {
    GT: pl.String,
    GQ: pl.Float64,
    PL_HOM_REF: pl.Float64,
    PL_HET: pl.Float64,
    PL_HOM_ALT: pl.Float64,
    OUTLIER_LOD: pl.Float64,
}
GENOTYPE_LIKELIHOOD_DTYPE: Final[pl.Struct] = pl.Struct(GENOTYPE_LIKELIHOOD_SCHEMA)

# Struct column the likelihood is computed into before being unnested
_GL: Final[str] = "_genotype_likelihood"

# Struct field carrying the per-observation ALT base error PROBABILITY into the
# kernel. Internal to the expression; it never reaches an output schema.
_ALT_BASE_ERROR: Final[str] = "_alt_base_error"


def _genotype_likelihood_batch(
    packed: pl.Series, eps: float, rho: float, theta_het: float
) -> pl.Series:
    """map_batches kernel: a (ref_reads, alt_reads, alt base error) struct Series
    in, a GENOTYPE_LIKELIHOOD_DTYPE struct Series out.

    The likelihood is a function of those three values alone, and they repeat
    enormously: with the downsample cap at 100 there are only a few thousand
    distinct combinations across tens of millions of cell-sites. So the special
    functions run once per distinct combination and the result is sent back. This
    is an exact deduplication, not an approximation, and it degrades gracefully
    when downsampling is off (worst case, every combination is distinct).

    Adding base quality to the key costs less than it looks: the instrument emits
    only four quality bins, so a group's mean error probability takes few values.
    Measured on the 10M run, the key goes from 3,347 distinct pairs to 9,080
    distinct triples over 656,611 rows, still a 72x reduction.
    """
    if packed.len() == 0:
        return pl.Series(_GL, [], dtype=GENOTYPE_LIKELIHOOD_DTYPE)

    fields = packed.struct.unnest()
    ref = fields[REF_READS].to_numpy().astype(np.int64)
    alt = fields[ALT_READS].to_numpy().astype(np.int64)
    # nan_to_num BEFORE the key, not just inside the likelihood: NaN never
    # compares equal to itself, so a single unmeasured row would otherwise defeat
    # the deduplication for every all-reference cell in the frame.
    err = np.nan_to_num(fields[_ALT_BASE_ERROR].to_numpy().astype(np.float64), nan=0.0)

    # Deduplicate on the whole (ref, alt, error) row. np.unique over the stacked
    # columns rather than one packed integer, since the error probability is a
    # float and cannot be folded into an integer key without losing exactness.
    _, first, inverse = np.unique(
        np.column_stack([ref, alt, err]), axis=0, return_index=True, return_inverse=True
    )
    inverse = inverse.reshape(-1)

    pl_homref, pl_het, pl_homalt, gq, gt_index, outlier_lod = (
        genotype_likelihood_arrays(
            ref[first],
            alt[first],
            eps=eps,
            rho=rho,
            theta_het=theta_het,
            alt_base_error=err[first],
        )
    )
    return pl.DataFrame(
        {
            GT: np.asarray(GT_LABELS)[gt_index][inverse],
            GQ: gq[inverse],
            PL_HOM_REF: pl_homref[inverse],
            PL_HET: pl_het[inverse],
            PL_HOM_ALT: pl_homalt[inverse],
            OUTLIER_LOD: outlier_lod[inverse],
        },
        schema=GENOTYPE_LIKELIHOOD_SCHEMA,
    ).to_struct(_GL)


def genotype_likelihood_expr(
    ref_col: str,
    alt_col: str,
    eps: float = DEFAULT_ERROR_EPS,
    rho: float = DEFAULT_ERROR_RHO,
    theta_het: float = THETA_HET,
    bq_col: str | None = None,
):
    """The genotype-likelihood struct expression, to be unnested by the caller.

    Null counts are read as zero: an all-reference cell has no ALT row to take a
    count from, and "no ALT reads observed" is exactly what the likelihood should
    be told. That is different from the quality statistics, where a null means
    "not measured" (see soft_filter_expr).

    bq_col names the per-allele mean base quality of the ALT-supporting reads and
    turns on the quality-aware hom-ref arm. The Phred is converted to an
    error probability HERE rather than in the kernel, so that the "not measured"
    null becomes an explicit 0.0 probability, which the floor then overrides. That
    keeps the two meanings of null apart: no ALT reads to measure gives 0.0 and
    falls back to the panel eps, which is not the same as a measured error of 0.
    Omit it (indels, which have no anchor base quality) for the flat-eps model.
    """
    alt_base_error = (
        pl.lit(0.0)
        if bq_col is None
        else (pl.lit(10.0) ** (-pl.col(bq_col) / 10.0)).fill_null(0.0)
    )
    return (
        pl.struct(
            pl.col(ref_col).fill_null(0).cast(pl.Int64).alias(REF_READS),
            pl.col(alt_col).fill_null(0).cast(pl.Int64).alias(ALT_READS),
            alt_base_error.cast(pl.Float64).alias(_ALT_BASE_ERROR),
        )
        .map_batches(
            lambda packed: _genotype_likelihood_batch(packed, eps, rho, theta_het),
            return_dtype=GENOTYPE_LIKELIHOOD_DTYPE,
            is_elementwise=True,
        )
        .alias(_GL)
    )


def genotype_filter_expr(gq_col: str, min_gq: float):
    """The genotype-row FILTER: genotype quality only, so PASS or lowGQ.

    Whether a genotype call is trustworthy is a question about how strongly the
    reads discriminate between the three genotypes, which is what GQ measures.
    Depth is subsumed rather than dropped: a cell with 3 reads mechanically
    cannot exceed GQ ~14, so it fails here without depth being named. Both DP
    and GQ stay reported, and the per-allele lowDP / lowVAF labels are
    unaffected on raw_variants (soft_filter_expr).

    A null GQ cannot arise from genotype_likelihood_expr (n = 0 scores GQ 0), so
    the fill_null here is defensive against an unjoined row only.
    """
    return (
        pl.when(pl.col(gq_col).fill_null(0.0) < min_gq)
        .then(pl.lit(FILTER_LOW_GQ))
        .otherwise(pl.lit(FILTER_PASS))
    )


def genotype_label_exprs(ref_col: str):
    """The (zygosity, variant_name) expression pair shared by filter_snps and
    filter_indels; they differ only in which column holds the reference allele.

    Both derive from GT, which is the argmin of the PLs. The VAF band that
    that used to decide this is deleted. Rows whose genotype FILTER is not PASS
    have too little evidence to call either way and get a null zygosity plus the
    lowqual name, which is the same contract as before so only the reason for
    being lowqual changed, from "too few reads" to "the reads do not separate
    the genotypes".

    variant_name keeps the historic strings so existing functions and calls (the h5ad obs
    layer) keep working: homo{ref}>{alt}, homo{ref}>{ref} for the reference
    state, and hete.
    """
    called = pl.col(FILTER) == FILTER_PASS
    hom_alt = pl.col(GT) == GT_HOM_ALT
    hom_ref = pl.col(GT) == GT_HOM_REF

    zygosity = (
        pl.when(~called)
        .then(pl.lit(None, dtype=pl.String))
        .when(hom_alt)
        .then(pl.lit(ZYG_HOM_ALT))
        .when(hom_ref)
        .then(pl.lit(ZYG_HOM_REF))
        .otherwise(pl.lit(ZYG_HET))
        .alias(ZYGOSITY)
    )
    variant_name = (
        pl.when(~called)
        .then(pl.lit(LOWQUAL_NAME))
        .when(hom_alt)
        .then(
            pl.concat_str(pl.lit("homo"), pl.col(ref_col), pl.lit(">"), pl.col(ALT_OUT))
        )
        .when(hom_ref)
        .then(
            pl.concat_str(pl.lit("homo"), pl.col(ref_col), pl.lit(">"), pl.col(ref_col))
        )
        .otherwise(pl.lit(HET_NAME))
        .alias(VARIANT_NAME)
    )
    return zygosity, variant_name


# Per-site cell-support counts. Since the cell barcode is the
# molecular unit at aggregation time, so one line of support will be "how many
# cells vote for the ALT", not just the raw-read VAF one over-amplified cell inflates.
CELLS_SUPPORTING_ALT: Final[str] = "cells_supporting_alt"
CELLS_TOTAL_AT_SITE: Final[str] = "cells_total_at_site"

# Genomic-coordinate columns (PR6 / COMB-532). Public output names: chrom is
# the real chromosome and pos is the 1-based genomic position. They live beside
# feature (the amplicon contig) and pos_local_0based (amplicon-local).
GENOMIC_CHROM: Final[str] = "chrom"
GENOMIC_POS: Final[str] = "pos"
# Self-describing FASTA header fields (bin/make_self_describing_fasta.py):
# g0 is the lift anchor, snp_pos the designed target site (optional).
G0: Final[str] = "g0"
TARGET_GENOMIC_POS: Final[str] = "snp_pos"
# Internal-only: the per-amplicon lift anchor carried from the FASTA header.
_G0: Final[str] = "_g0"

# Indel per-cell attributes (count_variant_indel.py). The allele identity is
# the indel signature carried per read out of the pileup, not a single base:
#
# 1. indel_len:
#   pysam PileupRead.indel at the anchor column (>0 INS, <0 DEL,
#   0 reference; null when the anchor base itself is deleted by a
#   larger event, so the read spans the locus but supports neither
#   the reference nor the called indel).
# 2. inserted_seq:
#   the inserted bases for an INS (query_sequence after the
#   anchor), "" for DEL/reference.
#
# INS matches the called allele on (indel_len AND inserted_seq)
# DEL matches on length alone (inserted_seq is "" on both the read and the called side, so the
# single join key (indel_len, inserted_seq) implements both rules).
INDEL_LEN: Final[str] = "indel_len"
INSERTED_SEQ: Final[str] = "inserted_seq"

# Self-describing per-cell indel feature key feature:pos:ref>alt
INDEL_FEATURE: Final[str] = "indel_feature"


# Output schema definition for the pileups
# Polars dtypes (UInt32 vs Int64 vs Int32) are explicitly set
@dataclass(frozen=True)
class ColumnSpec:
    """One output column: its name and its polars dtype."""

    name: str
    dtype: pl.DataType


def _empty_frame(schema: tuple[ColumnSpec, ...]) -> pl.DataFrame:
    """Zero-row DataFrame with exactly the columns/dtypes of schema."""
    return pl.DataFrame(schema={col.name: col.dtype for col in schema})


RAW_VARIANTS_EMPTY_SCHEMA: Final[tuple[ColumnSpec, ...]] = (
    ColumnSpec(BARCODE, pl.String),
    ColumnSpec(FEATURE, pl.String),
    ColumnSpec(POS_LOCAL, pl.Int64),
    ColumnSpec(GENOMIC_CHROM, pl.String),
    ColumnSpec(GENOMIC_POS, pl.Int32),
    ColumnSpec(REF_OUT, pl.String),
    ColumnSpec(ALT_OUT, pl.String),
    # UInt32 to match the live output. Otherwise an
    # Int64 declaration here makes the empty and
    # populated files fail to concat.
    ColumnSpec(REF_READS, pl.UInt32),
    ColumnSpec(ALT_READS, pl.UInt32),
    ColumnSpec(TOTAL_READS, pl.UInt32),
    ColumnSpec(VAF, pl.Float64),
    ColumnSpec(MEAN_BQ, pl.Float64),
    ColumnSpec(MEAN_MQ, pl.Float64),
    ColumnSpec(FILTER, pl.String),
    ColumnSpec(IS_NO_CALL, pl.Boolean),
    ColumnSpec(CALLERS, pl.List(pl.String)),
    ColumnSpec(N_CALLERS, pl.Int32),
    ColumnSpec(CELLS_SUPPORTING_ALT, pl.Int32),
    ColumnSpec(CELLS_TOTAL_AT_SITE, pl.Int32),
    ColumnSpec(SNP_FEATURE, pl.String),
)

RAW_INDELS_EMPTY_SCHEMA: Final[tuple[ColumnSpec, ...]] = (
    ColumnSpec(BARCODE, pl.String),
    ColumnSpec(FEATURE, pl.String),
    ColumnSpec(POS_LOCAL, pl.Int64),
    ColumnSpec(GENOMIC_CHROM, pl.String),
    ColumnSpec(GENOMIC_POS, pl.Int32),
    ColumnSpec(REF_OUT, pl.String),
    ColumnSpec(ALT_OUT, pl.String),
    ColumnSpec(REF_READS, pl.UInt32),
    ColumnSpec(ALT_READS, pl.UInt32),
    ColumnSpec(TOTAL_READS, pl.UInt32),
    ColumnSpec(VAF, pl.Float64),
    ColumnSpec(MEAN_MQ, pl.Float64),
    ColumnSpec(FILTER, pl.String),
    ColumnSpec(CALLERS, pl.List(pl.String)),
    ColumnSpec(N_CALLERS, pl.Int32),
    ColumnSpec(CELLS_SUPPORTING_ALT, pl.Int32),
    ColumnSpec(CELLS_TOTAL_AT_SITE, pl.Int32),
    ColumnSpec(INDEL_FEATURE, pl.String),
)


def pileup_truncated(bam, contig, start, stop, max_depth):
    """
    Obtain Pysam columns only at selected region.

    Base-quality filtering is intentionally disabled here (min_base_quality=0):
    this pileup is an attribution step that maps already-called variants back
    to cells, and the per-base quality filter belongs to the variant callers
    upstream (e.g., freebayes -q, bcftools mpileup -Q, GATK HC
    --base-quality-score-threshold, VarDict -q, etc). Re-applying a base-quality
    cutoff here would only discard reads the caller already accepted.
    """
    _, rtid, rstart, rstop = bam.parse_region(contig, start, stop)
    yield from IteratorColumnRegion(
        bam,
        tid=rtid,
        start=rstart,
        stop=rstop,
        truncate=True,
        max_depth=max_depth,
        min_base_quality=0,
    )


# TODO: get to actually only pileup on one position
# TODO: do a pileup by cell
def get_pileup_results(
    feature: str,
    start: int,
    stop: int,
    samfile_handle: pysam.AlignmentFile,
    reffa_handle: pysam.FastaFile,
    max_depth: int,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    # Iterate over pileup. Per read we carry the base quality at the anchor and the
    # read's mapping quality so the per-cell layer can *report* them (never gate on
    # them -- base-quality filtering stays with the callers upstream).

    results: list = []
    for pileupcolumn in pileup_truncated(
        samfile_handle,
        feature,
        start=start,
        stop=stop,
        max_depth=max_depth,
    ):
        pos = pileupcolumn.reference_pos  # 0-based
        if pos != start:
            continue
        chrom = pileupcolumn.reference_name
        # Get reference base (1-base to 1-base inclusive fetch)
        ref_base = reffa_handle.fetch(chrom, pos, pos + 1)  # returns a string of 1 base

        # Pooled site gate: total reads spanning the site, summed over every cell.
        if pileupcolumn.nsegments < min_site_pooled_coverage:
            continue
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                aln = pileupread.alignment
                read_base = aln.query_sequence[pileupread.query_position]
                # A read with no base qualities (QUAL "*") has query_qualities
                # None; report a null rather than raising, the mean just skips it.
                quals = aln.query_qualities
                base_qual = (
                    quals[pileupread.query_position] if quals is not None else None
                )
                map_qual = aln.mapping_quality
                results.append(
                    [aln.query_name, pos, ref_base, read_base, base_qual, map_qual]
                )
        if len(results) == 0:
            return None
        return pl.DataFrame(
            results,
            schema={
                READ_NAME: pl.String,
                "pos": pl.Int64,
                REF_BASE: pl.String,
                READ_BASE: pl.String,
                BASE_QUAL: pl.Int64,
                MAP_QUAL: pl.Int64,
            },
            orient="row",
        ).with_columns(pl.lit(feature).alias(FEATURE))
    return None


def build_read_cell_lut(
    barcode_mapping_path: str,
    filtered_amplicon_reads_path: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Builds a lookup table that resolves read to cell and requires them to map 1:1

    Returns (lut, cells):
      - lut: (read_name, cell_index), one row per barcoded read that is also in
        the amplicon allowlist.
      - cells: (barcode, cell_index), the reverse map, joined back onto the
        collapsed per-cell frame once the pileup loop is done.

    These are the same two inner joins the scripts used to run at read level
    after piling up every site: '.join(filtered_amplicon_reads).join(
    barcode_mapping)'
    """
    lut = (
        pl.scan_parquet(barcode_mapping_path)
        .select(READ_NAME, BARCODE)
        .join(
            pl.scan_parquet(filtered_amplicon_reads_path).select(READ_NAME),
            on=READ_NAME,
            how="inner",
        )
        .collect()
    )
    n_reads = lut[READ_NAME].n_unique()
    if lut.height != n_reads:
        raise ValueError(
            "FATAL ERROR:"
            f"read to cell map is not 1:1: {lut.height} rows for {n_reads} "
            f"reads, from {barcode_mapping_path} joined to "
            f"{filtered_amplicon_reads_path}. Each repeat would be counted "
            "again at every site the read spans; deduplicate whichever input "
            "repeats, not this table."
        )
    cells = lut.select(BARCODE).unique().sort(BARCODE).with_row_index(CELL_INDEX)
    return lut.join(cells, on=BARCODE).select(READ_NAME, CELL_INDEX), cells


def aggregate_snv_sites(
    contig: str,
    sites: pl.DataFrame,
    samfile_handle: pysam.AlignmentFile,
    ref_handle: pysam.FastaFile,
    lut: pl.DataFrame,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    """Pile up every selected column of one contig collapsed per cell.

    'sites' holds the 0-based columns to pile up on 'contig'. The collapse runs
    per site, inside the loop, so the widest thing in memory is one site's reads
    rather than every (read, site) pair of the sample. What leaves this function
    is already (cell x allele). Reads that belong to no called cell are dropped
    by the inner join on the lut. None when no site on this contig produced a
    single such read.
    """
    max_depth = samfile_handle.mapped
    per_site: list[pl.DataFrame] = []
    for start in sites[POS]:
        reads = get_pileup_results(
            contig,
            start,
            stop=start + 1,
            samfile_handle=samfile_handle,
            reffa_handle=ref_handle,
            max_depth=max_depth,
            min_site_pooled_coverage=min_site_pooled_coverage,
        )
        if reads is None or reads.shape[0] == 0:
            continue
        counts = format_raw_snps(
            reads.lazy().join(lut.lazy(), on=READ_NAME).rename({"pos": POS}),
            cell_col=CELL_INDEX,
        ).collect()
        if counts.shape[0] == 0:
            continue
        per_site.append(counts)
    if not per_site:
        return None
    return pl.concat(per_site)


def per_cell_snv_counts(
    selected_features: pl.DataFrame,
    samfile_handle: pysam.AlignmentFile,
    ref_handle: pysam.FastaFile,
    lut: pl.DataFrame,
    cells: pl.DataFrame,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    """Per-cell per-allele read counts over every selected SNV site.

    One aggregate_snv_sites call per contig, then the barcode is attached once,
    on the collapsed dataframe. Returns None when no site produced any read from a
    called cell, so the caller will short-circuit to the empty-output path,
    mirroring per_cell_indel_counts.

    Grouping the sites by contig is free O(reads x sites) and this does not change it.
    """
    per_contig: list[pl.DataFrame] = []
    for (contig,), sites in selected_features.group_by(FEATURE, maintain_order=True):
        counts = aggregate_snv_sites(
            contig,
            sites,
            samfile_handle=samfile_handle,
            ref_handle=ref_handle,
            lut=lut,
            min_site_pooled_coverage=min_site_pooled_coverage,
        )
        if counts is not None:
            per_contig.append(counts)
    if not per_contig:
        return None
    return (
        pl.concat(per_contig)
        .join(cells, on=CELL_INDEX)
        .select(
            BARCODE,
            POS,
            REF_BASE,
            FEATURE,
            READS,
            MEAN_BQ,
            MEAN_MQ,
            TOTAL_READS,
            BASES,
            SNP_FEATURE,
        )
    )


def format_raw_snps(
    raw_barcoded_snps: pl.LazyFrame, cell_col: str = BARCODE
) -> pl.LazyFrame:
    """Collapse per-read SNV observations to per-cell per-allele counts.

    One row per (cell, pos, feature, allele) with the read count, the quality
    aggregates of the supporting reads, the cell-site total broadcast over every
    allele, and the ref>alt / snp_feature labels.

    cell_col names the column that identifies the cell. The pileup loop
    aggregates on the dense CELL_INDEX and attaches the barcode afterwards
    (see per_cell_snv_counts); a data frame that already carries barcodes
    takes the default.
    """

    df = (
        raw_barcoded_snps.group_by(cell_col, POS, READ_BASE, REF_BASE, FEATURE)
        .agg(
            pl.len().alias(READS),
            # per-allele base / mapping quality of the supporting reads. Base
            # quality is averaged in probability space and returned to Phred
            # (mean_bq_expr) so mapping quality stays a plain mean, since lowMQ is
            # read about a group of alignments rather than fed to a likelihood.
            mean_bq_expr(BASE_QUAL),
            pl.col(MAP_QUAL).mean().alias(MEAN_MQ),
        )
        .with_columns(pl.sum(READS).over(cell_col, POS, FEATURE).alias(TOTAL_READS))
    )
    raw_snp = (
        df.with_columns(
            bases=pl.concat_str(pl.col(REF_BASE), pl.col(READ_BASE), separator=">")
        )
        .with_columns(pl.format("{}_{}_{}", FEATURE, POS, READ_BASE).alias(SNP_FEATURE))
        .drop(READ_BASE)
    )
    return raw_snp


def _optional_header_int(value: str | None) -> int | None:
    """An optional integer FASTA-header field, where unparseable reads as absent.

    Headers are hand-editable, so `snp_pos=`, `snp_pos=NA` and a stray space
    after the `=` all occur.

    Out-of-range values are rejected too, since the column is Int32.

    TODO [COMB-585]: the header specification allows a pipe-separated list, `snp_pos=250|260`, and
    the FASTA accepts one. Here `int()` raises on it and the field reads as absent,
    so such an amplicon silently gets a null target position and its on-target
    classification is wrong with no error anywhere. Either take every value or narrow
    the specification to a single integer; today the two disagree.
    """
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if -(2**31) <= parsed < 2**31 else None


def parse_amplicon_headers(fasta_path: str) -> pl.DataFrame:
    """Parse the self-describing amplicon FASTA headers into one row per contig

    The read is intentionally dependency-free: pysam.FastaFile generates record
    IDs only (it truncates at the first whitespace), not the description where the
    metadata lives, so the > lines are parsed directly. The locked feature is
    the header format (see bin/make_self_describing_fasta.py), not this parser.

    A plain (non-enriched) FASTA will simply generate an empty table.

    TODO [COMB-586]: keys are matched case-sensitively here, while the FASTA check lowercases them
    before checking. A header written `Chrom=` or `G0=` therefore validates clean and is
    then silently skipped below, leaving null chrom and position with no error anywhere.
    Either read the headers through `lib.fasta.parse_fasta`, which lowercases, or drop the
    lowercasing there so the gate matches this parser.
    """
    feats: list[str] = []
    chroms: list[str] = []
    g0s: list[int] = []
    targets: list[int | None] = []
    with open(fasta_path, encoding="utf-8-sig") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            tokens = line[1:].split()
            fields = dict(t.split("=", 1) for t in tokens[1:] if "=" in t)
            if GENOMIC_CHROM not in fields or G0 not in fields:
                continue  # plain header without genomic metadata
            feats.append(tokens[0])
            chroms.append(fields[GENOMIC_CHROM])
            g0s.append(int(fields[G0]))
            targets.append(_optional_header_int(fields.get(TARGET_GENOMIC_POS)))
    return pl.DataFrame(
        {
            FEATURE: feats,
            GENOMIC_CHROM: chroms,
            G0: g0s,
            TARGET_GENOMIC_POS: targets,
        },
        schema={
            FEATURE: pl.String,
            GENOMIC_CHROM: pl.String,
            G0: pl.Int32,
            TARGET_GENOMIC_POS: pl.Int32,
        },
    )


def parse_amplicon_genomic_map(fasta_path: str) -> pl.DataFrame:
    """The amplicon to genome lift table used by add_genomic_coords: feature,
    genomic chrom and the lift anchor, the latter under the internal _g0 name so
    it cannot collide with an output column.

    A plain (non-enriched) FASTA will simply generate an empty table. The lift then adds
    null chrom/pos, so the pipeline finishes gracefully rather than failing.
    """
    return parse_amplicon_headers(fasta_path).select(
        FEATURE, GENOMIC_CHROM, pl.col(G0).alias(_G0)
    )


def add_genomic_coords(frame: pl.LazyFrame, genomic_map: pl.DataFrame) -> pl.LazyFrame:
    """Left-join the amplicon to genome lift onto a per-cell framework and then add the
    1-based genomic chrom/pos.

    pos = g0 + pos_local_0based: g0 is the 1-based genomic coordinate aligned
    to amplicon-local 0-based position 0, so there is no weird +1. Features absent
    from the map get null chrom/pos.
    """
    return (
        frame.join(genomic_map.lazy(), on=FEATURE, how="left")
        .with_columns(
            (pl.col(_G0) + pl.col(POS_LOCAL)).cast(pl.Int32).alias(GENOMIC_POS)
        )
        .drop(_G0)
    )


def per_site_cell_support(
    raw_snp: pl.LazyFrame,
    flag_min_dp: int,
    flag_min_vaf: float,
) -> pl.LazyFrame:
    """Per-site cell-support counts for the pseudobulk calls

    Two counts per called allele is calculated grouped on (feature, pos_local_0based, ref, alt)
    so they join exactly like callers/n_callers:

      - cells_total_at_site: cells with callable coverage at the site
        (total_reads >= flag_min_dp); the denominator across every allele of the site.
      - cells_supporting_alt: of those, the cells whose per-cell ALT fraction
        passes flag_min_vaf for this allele (a real ALT, not the ref base and not
        the N). These are the same flag thresholds count_variant labels a per-cell
        call PASS with, so the count equals the cells that would PASS for the allele.

    The soft FILTER never deletes rows, but the support counts still need a
    definition of "callable" and "supporting"; the flag thresholds provide it.
    Ref and N rows will have no row and so join to null downstream, again it is exactly
    like callers/n_callers for alleles that no caller reported them.
    """
    eligible = raw_snp.filter(pl.col(TOTAL_READS) >= flag_min_dp)

    cells_total = eligible.group_by(FEATURE, POS).agg(
        pl.col(BARCODE).n_unique().cast(pl.Int32).alias(CELLS_TOTAL_AT_SITE)
    )

    cells_alt = (
        eligible.with_columns(pl.col(BASES).str.split(">").list.last().alias(ALT_OUT))
        .filter(
            (pl.col(ALT_OUT) != pl.col(REF_BASE))
            & (pl.col(ALT_OUT) != "N")
            & ((pl.col(READS) / pl.col(TOTAL_READS)) >= flag_min_vaf)
        )
        .group_by(FEATURE, POS, pl.col(REF_BASE).alias(REF_OUT), ALT_OUT)
        .agg(pl.col(BARCODE).n_unique().cast(pl.Int32).alias(CELLS_SUPPORTING_ALT))
    )

    return (
        cells_alt.join(cells_total, on=[FEATURE, POS], how="left")
        .rename({POS: POS_LOCAL})
        .select(
            FEATURE,
            POS_LOCAL,
            REF_OUT,
            ALT_OUT,
            CELLS_SUPPORTING_ALT,
            CELLS_TOTAL_AT_SITE,
        )
    )


def format_raw_snps_wide(
    raw_snp: pl.LazyFrame,
    callers_support: pl.LazyFrame,
    cell_support: pl.LazyFrame,
    flag_min_dp: int,
    flag_min_vaf: float,
    flag_min_bq: float = DEFAULT_FLAG_MIN_BQ,
    flag_min_mq: float = DEFAULT_FLAG_MIN_MQ,
    genomic_map: pl.DataFrame | None = None,
) -> pl.LazyFrame:
    """Enrich the long per-allele raw SNP frame into the downstream-friendly
    raw_variants.parquet schema.

    Adds:
      - explicit ref / alt (alt is the per-row allele; ref is the site base)
      - wide ref_reads (the site's reference-supporting count, broadcast),
        alt_reads (this allele's count) and vaf = alt_reads / total_reads
      - mean_bq / mean_mq: per-allele mean base / mapping quality (reported only)
      - filter: the soft per-allele label (PASS, or the semicolon-joined
        lowDP / lowVAF / lowBQ / lowMQ reasons) for this allele row; nothing is
        dropped, the label just annotates it
      - is_no_call for the N allele rows (denominator left unchanged)
      - callers (sorted list of callers backing the allele) and its length
        n_callers, joined per (feature, pos, ref, alt); both are null for
        alleles no caller reported (ref rows, N, sub-threshold noise)
      - cells_supporting_alt / cells_total_at_site joined on the
        same (feature, pos, ref, alt) grouping from per_site_cell_support
      - the local coordinate renamed to the self-documenting pos_local_0based
      - genomic chrom/pos (1-based) when genomic_map is supplied; placed
        right after pos_local_0based. Omitted when no map is given (e.g. a
        plain FASTA ), keeping the call backward
        compatible.
    """
    wide = (
        raw_snp.with_columns(
            pl.col(REF_BASE).alias(REF_OUT),
            pl.col(BASES).str.split(">").list.last().alias(ALT_OUT),
            pl.col(READS).alias(ALT_READS),
        )
        .with_columns(
            (pl.col(ALT_OUT) == "N").alias(IS_NO_CALL),
            pl.when(pl.col(ALT_OUT) == pl.col(REF_OUT))
            .then(pl.col(ALT_READS))
            .otherwise(0)
            .sum()
            .over(BARCODE, FEATURE, POS)
            .alias(REF_READS),
        )
        .with_columns((pl.col(ALT_READS) / pl.col(TOTAL_READS)).alias(VAF))
        .with_columns(
            soft_filter_expr(
                TOTAL_READS,
                VAF,
                flag_min_dp,
                flag_min_vaf,
                bq_col=MEAN_BQ,
                mq_col=MEAN_MQ,
                flag_min_bq=flag_min_bq,
                flag_min_mq=flag_min_mq,
            ).alias(FILTER)
        )
        .rename({POS: POS_LOCAL})
        .join(callers_support, on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT], how="left")
        .join(cell_support, on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT], how="left")
    )

    leading = [BARCODE, FEATURE, POS_LOCAL]
    if genomic_map is not None:
        wide = add_genomic_coords(wide, genomic_map)
        leading += [GENOMIC_CHROM, GENOMIC_POS]

    return wide.select(
        *leading,
        REF_OUT,
        ALT_OUT,
        REF_READS,
        ALT_READS,
        TOTAL_READS,
        VAF,
        MEAN_BQ,
        MEAN_MQ,
        FILTER,
        IS_NO_CALL,
        CALLERS,
        N_CALLERS,
        CELLS_SUPPORTING_ALT,
        CELLS_TOTAL_AT_SITE,
        SNP_FEATURE,
    )


def empty_raw_variants() -> pl.DataFrame:
    """Zero-row raw_variants.parquet written by count_variant.py when a sample has
    no called positions. Schema is RAW_VARIANTS_EMPTY_SCHEMA.
    """
    return _empty_frame(RAW_VARIANTS_EMPTY_SCHEMA)


def filter_snps(
    raw_snp: pl.LazyFrame,
    min_gq: float,
    eps: float = DEFAULT_ERROR_EPS,
    rho: float = DEFAULT_ERROR_RHO,
) -> pl.LazyFrame:
    """Per-cell genotype for each covered cell-site soft-filtered only, nothing is dropped.

    Reduces the long per-allele frame to one row per (barcode, feature, pos_local):
    the genotype ALT is the highest-VAF non-reference, non-N allele in the cell
    (ties broken by base), and every covered cell-site is kept, including
    all-reference cells (which come out as a called homref).

    The genotype itself comes from the beta-binomial likelihood over
    (ref_reads, alt_reads, mean_bq) - please see the genotype_likelihood_arrays.
    This is so the row carries GT, GQ, the three PLs and outlier_lod alongside its
    depth, VAF and mean base/mapping quality. Nothing is dropped for quality: a
    cell whose reads do not separate the genotypes is labelled (FILTER = lowGQ,
    variant_name = lowqual) rather than deleted.

    The chosen ALT's own base quality raises the hom-ref arm, so a well-covered reference
    cell holding one junk read stays a called homref instead of falling to lowqual. It is
    the same statistic the lowBQ label reads, used for a different question, and it is taken
    from the chosen ALT allele only. See the notes in the genotype_likelihood_arrays function
    above.
    """
    alleles = raw_snp.with_columns(
        pl.col(BASES).str.split(">").list.last().alias(ALT_OUT),
        (pl.col(READS) / pl.col(TOTAL_READS)).alias(VAF),
    )

    # Every covered cell-site is the denominator; keep all-ref / all-N cells too.
    # ref_reads is the site's reference-supporting count in this cell, which the
    # likelihood needs as the other half of n (third alleles and N belong to
    # neither side and are intentionally left out of it).
    sites = (
        alleles.with_columns(
            pl.when(pl.col(ALT_OUT) == pl.col(REF_BASE))
            .then(pl.col(READS))
            .otherwise(0)
            .sum()
            .over(BARCODE, FEATURE, POS)
            .alias(REF_READS)
        )
        .select(BARCODE, POS, FEATURE, REF_BASE, TOTAL_READS, REF_READS)
        .unique()
    )

    # Genotype ALT = highest-VAF real allele (non-ref, non-N) in the cell; ties
    # broken by base so the pick is deterministic regardless of scan/join order.
    chosen = (
        alleles.filter((pl.col(ALT_OUT) != pl.col(REF_BASE)) & (pl.col(ALT_OUT) != "N"))
        .sort([VAF, ALT_OUT], descending=[True, False])
        .group_by([BARCODE, POS, FEATURE], maintain_order=True)
        .agg(
            pl.col(ALT_OUT).first(),
            pl.col(READS).first().alias(ALT_READS),
            pl.col(VAF).first(),
            pl.col(MEAN_BQ).first(),
            pl.col(MEAN_MQ).first(),
        )
    )

    return (
        sites.join(chosen, on=[BARCODE, POS, FEATURE], how="left")
        .with_columns(
            pl.col(ALT_READS).fill_null(0),
            pl.col(VAF).fill_null(0.0),
        )
        .with_columns(
            genotype_likelihood_expr(
                REF_READS, ALT_READS, eps=eps, rho=rho, bq_col=MEAN_BQ
            )
        )
        .unnest(_GL)
        .with_columns(genotype_filter_expr(GQ, min_gq).alias(FILTER))
        .with_columns(*genotype_label_exprs(REF_BASE))
        .rename({POS: POS_LOCAL, REF_BASE: REF_OUT})
        .select(
            BARCODE,
            POS_LOCAL,
            FEATURE,
            REF_OUT,
            ALT_OUT,
            TOTAL_READS,
            REF_READS,
            ALT_READS,
            VAF,
            MEAN_BQ,
            MEAN_MQ,
            GT,
            GQ,
            PL_HOM_REF,
            PL_HET,
            PL_HOM_ALT,
            OUTLIER_LOD,
            ZYGOSITY,
            VARIANT_NAME,
            FILTER,
        )
    )


#############################
# Indel per-cell attribution (count_variant_indel.py).
#
# The SNV helpers above read a single base per read and everything is on a
# REF>ALT base string; get_pileup_results deliberately skips is_del
# reads. An indel can _never_ be represented that way, so the functions below take
# a parallel path: they read pysam.PileupRead.indel at the anchor column
# on the indel signature (indel_len, inserted_seq). They share the
# coordinate lift (parse_amplicon_genomic_map / add_genomic_coords) and the
# homo/hete thresholds with the SNV path; nothing here mutates the SNV output.
#############################


def get_indel_pileup_results(
    feature: str,
    start: int,
    stop: int,
    samfile_handle: pysam.AlignmentFile,
    max_depth: int,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    """Pileup the single anchor column (0-based start) and generate one row per
    spanning read: (read_name, pos, indel_len, inserted_seq, map_qual).

    Unlike get_pileup_results this will is_del or no-query reads so the site
    stays "reads spanning the anchor". Those reads carry a null
    indel_len so they match neither the reference (indel_len == 0) nor any
    called indel allele. For an insertion (indel_len > 0) the inserted bases
    are read from the query sequence immediately after the anchor, deletions and
    reference reads carry inserted_seq == "". The read's mapping quality is carried
    so the per-cell layer can report it (indels have no single anchor base quality).
    """
    for pileupcolumn in pileup_truncated(
        samfile_handle,
        feature,
        start=start,
        stop=stop,
        max_depth=max_depth,
    ):
        pos = pileupcolumn.reference_pos  # 0-based!
        if pos != start:
            continue
        if pileupcolumn.nsegments < min_site_pooled_coverage:
            return None
        results: list = []
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_refskip:
                continue
            aln = pileupread.alignment
            map_qual = aln.mapping_quality
            query_position = pileupread.query_position
            if pileupread.is_del or query_position is None:
                # NOTE: an edge case where the anchor base is deleted by a larger event
                # This large even can span the locus but supports neither reference nor the called indel.
                results.append([aln.query_name, pos, None, "", map_qual])
                continue
            indel_len = pileupread.indel
            if indel_len > 0:
                inserted_seq = aln.query_sequence[
                    query_position + 1 : query_position + 1 + indel_len
                ]
            else:
                inserted_seq = ""
            results.append([aln.query_name, pos, indel_len, inserted_seq, map_qual])
        if len(results) == 0:
            return None

        # Add explicit schema
        return pl.DataFrame(
            results,
            schema={
                READ_NAME: pl.String,
                "pos": pl.Int64,
                INDEL_LEN: pl.Int64,
                INSERTED_SEQ: pl.String,
                MAP_QUAL: pl.Int64,
            },
            orient="row",
        ).with_columns(pl.lit(feature).alias(FEATURE))
    return None


def aggregate_indel_sites(
    contig: str,
    sites: pl.DataFrame,
    samfile_handle: pysam.AlignmentFile,
    lut: pl.DataFrame,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    """Pile up every selected anchor column of one contig, collapsed per cell.

    The indel mirror of aggregate_snv_sites: the collapse runs per anchor,
    inside the loop, so the frame that leaves this function is (cell x
    signature) and never (read x site). None when no anchor on this contig
    produced a spanning read from a called cell.
    """
    max_depth = samfile_handle.mapped
    per_site: list[pl.DataFrame] = []
    for start in sites[POS]:
        reads = get_indel_pileup_results(
            contig,
            start,
            stop=start + 1,
            samfile_handle=samfile_handle,
            max_depth=max_depth,
            min_site_pooled_coverage=min_site_pooled_coverage,
        )
        if reads is None or reads.shape[0] == 0:
            continue
        counts = format_raw_indels(
            reads.lazy().join(lut.lazy(), on=READ_NAME).rename({"pos": POS_LOCAL}),
            cell_col=CELL_INDEX,
        ).collect()
        if counts.shape[0] == 0:
            continue
        per_site.append(counts)
    if not per_site:
        return None
    return pl.concat(per_site)


def per_cell_indel_counts(
    selected_features: pl.DataFrame,
    samfile_handle: pysam.AlignmentFile,
    lut: pl.DataFrame,
    cells: pl.DataFrame,
    min_site_pooled_coverage: int,
) -> pl.DataFrame | None:
    """Per-cell per-signature read counts over every selected indel anchor.

    One aggregate_indel_sites call per contig, then the barcode is attached
    once, on the collapsed frame. Returns None when no anchor produced any
    spanning read from a called cell, so the caller can short-circuit to the
    empty-output path.
    """
    per_contig: list[pl.DataFrame] = []
    for (contig,), sites in selected_features.group_by(FEATURE, maintain_order=True):
        counts = aggregate_indel_sites(
            contig,
            sites,
            samfile_handle=samfile_handle,
            lut=lut,
            min_site_pooled_coverage=min_site_pooled_coverage,
        )
        if counts is not None:
            per_contig.append(counts)
    if not per_contig:
        return None
    return (
        pl.concat(per_contig)
        .join(cells, on=CELL_INDEX)
        .select(
            BARCODE,
            POS_LOCAL,
            FEATURE,
            INDEL_LEN,
            INSERTED_SEQ,
            READS,
            MEAN_MQ,
            TOTAL_READS,
        )
    )


def format_raw_indels(
    raw_barcoded_indels: pl.LazyFrame, cell_col: str = BARCODE
) -> pl.LazyFrame:
    """Collapse per-read indel observations to per-cell per-signature counts.

    One row per (cell, pos_local, feature, indel_len, inserted_seq) with the
    read count, the mean mapping quality of the supporting reads, and the site
    total (every spanning read in the cell, broadcast over the signature).
    Mirrors format_raw_snps above, including cell_col, but the grouping key is
    the indel signature.
    """
    return (
        raw_barcoded_indels.group_by(
            cell_col, POS_LOCAL, FEATURE, INDEL_LEN, INSERTED_SEQ
        )
        .agg(
            pl.len().alias(READS),
            pl.col(MAP_QUAL).mean().alias(MEAN_MQ),
        )
        .with_columns(
            pl.col(READS).sum().over(cell_col, POS_LOCAL, FEATURE).alias(TOTAL_READS)
        )
    )


def indel_called_alleles(calls: pl.LazyFrame) -> pl.LazyFrame:
    """Per called INS/DEL allele: the consensus support and the join signature.

    'calls' is the catalog parquet already filtered to 'type IN ('INS','DEL')'.
    Groups to one row per (feature, pos_local_0based, ref, alt) carrying the
    sorted callers list / n_callers, plus the (indel_len, inserted_seq)
    signature used to match reads:
      - INS: indel_len = len(alt) - len(ref) > 0, inserted_seq = alt[len(ref):]
      - DEL: indel_len = -(len(ref) - len(alt)) < 0, inserted_seq = ""
    POS is 1-based in the catalog, so -1 lifts it to the 0-based anchor column.
    """
    # len_bytes() is UInt32; cast to signed before subtracting so a deletion
    # (len_alt < len_ref) yields a negative indel_len instead of underflowing.
    len_ref = pl.col(REF_OUT).str.len_bytes().cast(pl.Int64)
    len_alt = pl.col(ALT_OUT).str.len_bytes().cast(pl.Int64)
    return (
        calls.group_by(
            pl.col("chrom").alias(FEATURE),
            (pl.col("pos").cast(pl.Int64) - 1).alias(POS_LOCAL),
            pl.col("ref").alias(REF_OUT),
            pl.col("alt").alias(ALT_OUT),
        )
        .agg(
            pl.col("caller").unique().sort().alias(CALLERS),
            pl.col("caller").n_unique().cast(pl.Int32).alias(N_CALLERS),
        )
        .with_columns(
            (len_alt - len_ref).cast(pl.Int64).alias(INDEL_LEN),
            pl.when(len_alt > len_ref)
            .then(pl.col(ALT_OUT).str.slice(len_ref))
            .otherwise(pl.lit(""))
            .alias(INSERTED_SEQ),
        )
    )


def per_site_indel_cell_support(
    observed: pl.LazyFrame,
    alt_long: pl.LazyFrame,
    flag_min_dp: int,
    flag_min_vaf: float,
) -> pl.LazyFrame:
    """Per-site cell-support counts for the indel calls, the indel analogue of
    per_site_cell_support:

      - cells_total_at_site: cells with callable coverage at the anchor
        (total_reads >= flag_min_dp), counted over every spanning read so
        it is the denominator shared by every allele of the site.
      - cells_supporting_alt: of those, the cells whose per-cell VAF for this
        called allele passes flag_min_vaf.

    Same flag thresholds count_variant_indel labels a per-cell call PASS with.
    Both join downstream on (feature, pos_local, ref, alt) exactly like
    callers/n_callers.
    """
    cells_total = (
        observed.filter(pl.col(TOTAL_READS) >= flag_min_dp)
        .group_by(FEATURE, POS_LOCAL)
        .agg(pl.col(BARCODE).n_unique().cast(pl.Int32).alias(CELLS_TOTAL_AT_SITE))
    )
    cells_alt = alt_long.group_by(FEATURE, POS_LOCAL, REF_OUT, ALT_OUT).agg(
        pl.col(BARCODE)
        .filter(
            (pl.col(ALT_READS) > 0)
            & (pl.col(TOTAL_READS) >= flag_min_dp)
            & (pl.col(VAF) >= flag_min_vaf)
        )
        .n_unique()
        .cast(pl.Int32)
        .alias(CELLS_SUPPORTING_ALT)
    )
    return cells_alt.join(cells_total, on=[FEATURE, POS_LOCAL], how="left").select(
        FEATURE,
        POS_LOCAL,
        REF_OUT,
        ALT_OUT,
        CELLS_SUPPORTING_ALT,
        CELLS_TOTAL_AT_SITE,
    )


def indel_anchor_coverage(observed: pl.LazyFrame) -> pl.LazyFrame:
    """Per (barcode, feature, anchor): anchor depth and reference-supporting reads.

    An anchor-deleted read has a null indel_len that no comparison matches, so
    it lands outside ref_reads and outside every allele's alt_reads: that is the
    DP > 0, REF + AD == 0 cell with nothing to genotype.
    """
    return observed.group_by(BARCODE, FEATURE, POS_LOCAL).agg(
        pl.col(TOTAL_READS).first(),
        pl.col(READS).filter(pl.col(INDEL_LEN) == 0).sum().alias(REF_READS),
    )


def indel_alt_long(
    observed: pl.LazyFrame,
    called_alleles: pl.LazyFrame,
) -> pl.LazyFrame:
    """One row per (cell spanning a called anchor, called allele at that anchor).

    Coverage drives the join, not the observed signature, so a covered
    non-carrier arrives at alt_reads=0 and filter_indels calls it hom-ref.

    Per EVENT, never per anchor: events share indel anchors in the VCF sense
    and a carrier of one is a non-carrier of the other, so summing total_reads
    over rows double counts a shared anchor.

    total_reads is the wider denominator (anchor-deleted reads included);
    ref_reads + alt_reads is the n the likelihood needs. mean_mq stays null for
    a non-carrier, which soft_filter_expr reads as unknown, not as a failure.
    """
    return (
        indel_anchor_coverage(observed)
        .join(called_alleles, on=[FEATURE, POS_LOCAL], how="inner")
        .join(
            observed.select(
                BARCODE, FEATURE, POS_LOCAL, INDEL_LEN, INSERTED_SEQ, READS, MEAN_MQ
            ),
            on=[BARCODE, FEATURE, POS_LOCAL, INDEL_LEN, INSERTED_SEQ],
            how="left",
            validate="m:1",
        )
        .with_columns(pl.col(READS).fill_null(0).alias(ALT_READS))
        .with_columns(
            (pl.col(ALT_READS) / pl.col(TOTAL_READS)).alias(VAF),
            pl.format("{}:{}:{}>{}", FEATURE, POS_LOCAL, REF_OUT, ALT_OUT).alias(
                INDEL_FEATURE
            ),
        )
    )


def format_raw_indels_wide(
    alt_long: pl.LazyFrame,
    cell_support: pl.LazyFrame,
    flag_min_dp: int,
    flag_min_vaf: float,
    flag_min_mq: float = DEFAULT_FLAG_MIN_MQ,
    genomic_map: pl.DataFrame | None = None,
) -> pl.LazyFrame:
    """An indel sibling of format_raw_snps_wide

    ref_reads (the cell's reference-supporting count at the anchor), mean_mq,
    callers/n_callers and the cell-support counts all ride along from alt_long /
    cell_support. A soft filter label (PASS, or the joined lowDP / lowVAF / lowMQ
    reasons) is added per allele row; nothing is dropped. lowBQ cannot appear
    here: an indel has no single anchor base quality, so no mean_bq is measured.
    Genomic chrom/pos are added when a genomic_map is supplied.
    """
    wide = alt_long.join(
        cell_support, on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT], how="left"
    ).with_columns(
        soft_filter_expr(
            TOTAL_READS,
            VAF,
            flag_min_dp,
            flag_min_vaf,
            mq_col=MEAN_MQ,
            flag_min_mq=flag_min_mq,
        ).alias(FILTER)
    )
    leading = [BARCODE, FEATURE, POS_LOCAL]
    if genomic_map is not None:
        wide = add_genomic_coords(wide, genomic_map)
        leading += [GENOMIC_CHROM, GENOMIC_POS]
    return wide.select(
        *leading,
        REF_OUT,
        ALT_OUT,
        REF_READS,
        ALT_READS,
        TOTAL_READS,
        VAF,
        MEAN_MQ,
        FILTER,
        CALLERS,
        N_CALLERS,
        CELLS_SUPPORTING_ALT,
        CELLS_TOTAL_AT_SITE,
        INDEL_FEATURE,
    )


def unattributed_indel_site_rows(
    called_alleles: pl.LazyFrame,
    alt_long: pl.LazyFrame,
    genomic_map: pl.DataFrame | None = None,
) -> pl.LazyFrame:
    """Site-level placeholder rows for called INS/DEL alleles no cell covers.

    Now that indel_alt_long is coverage-driven, "absent from alt_long" means the
    anchor was never piled up, which is the only class left with no per-cell
    trace. An allele that WAS covered and simply carried by nobody describes
    itself through its own rows, so it must not get a placeholder too: it would
    be the same allele reported twice, once at cells_supporting_alt 0 and once
    per cell.

    Here we anti-join those alleles back and create one zero-support placeholder each
    (barcode null, alt_reads/cell counts 0, reads/vaf null) so nothing silently
    disappears; the schema matches format_raw_indels_wide so the two concat
    cleanly. Un-pileupable alleles then surface at cells_supporting_alt = 0 (and at
    y=0 in the abundance plot) and become countable for the indels_attributed / indels_called
    diagnostic.
    """
    dropped = called_alleles.join(
        alt_long.select(FEATURE, POS_LOCAL, REF_OUT, ALT_OUT).unique(),
        on=[FEATURE, POS_LOCAL, REF_OUT, ALT_OUT],
        how="anti",
    ).with_columns(
        pl.lit(None, dtype=pl.String).alias(BARCODE),
        pl.lit(None, dtype=pl.UInt32).alias(REF_READS),
        pl.lit(0, dtype=pl.UInt32).alias(ALT_READS),
        pl.lit(None, dtype=pl.UInt32).alias(TOTAL_READS),
        pl.lit(None, dtype=pl.Float64).alias(VAF),
        pl.lit(None, dtype=pl.Float64).alias(MEAN_MQ),
        # No read attributed, so no per-cell call to label.
        pl.lit(None, dtype=pl.String).alias(FILTER),
        pl.lit(0, dtype=pl.Int32).alias(CELLS_SUPPORTING_ALT),
        pl.lit(0, dtype=pl.Int32).alias(CELLS_TOTAL_AT_SITE),
        pl.format("{}:{}:{}>{}", FEATURE, POS_LOCAL, REF_OUT, ALT_OUT).alias(
            INDEL_FEATURE
        ),
    )

    leading = [BARCODE, FEATURE, POS_LOCAL]
    if genomic_map is not None:
        dropped = add_genomic_coords(dropped, genomic_map)
        leading += [GENOMIC_CHROM, GENOMIC_POS]

    return dropped.select(
        *leading,
        REF_OUT,
        ALT_OUT,
        REF_READS,
        ALT_READS,
        TOTAL_READS,
        VAF,
        MEAN_MQ,
        FILTER,
        CALLERS,
        N_CALLERS,
        CELLS_SUPPORTING_ALT,
        CELLS_TOTAL_AT_SITE,
        INDEL_FEATURE,
    )


def empty_raw_indels() -> pl.DataFrame:
    """Zero-row raw_variants_indel.parquet written by count_variant_indel.py when
    a sample has no called indels (or no spanning reads), mirroring
    empty_raw_variants that we have for the SNV path.
    """
    return _empty_frame(RAW_INDELS_EMPTY_SCHEMA)


def filter_indels(
    alt_long: pl.LazyFrame,
    min_gq: float,
    eps: float = DEFAULT_ERROR_EPS,
    rho: float = DEFAULT_ERROR_RHO,
) -> pl.LazyFrame:
    """Per-cell genotype for each called indel allele -- soft-filtered, nothing
    dropped; the indel analogue of filter_snps.

    Every (barcode, called allele) row is kept and genotyped from the same
    beta-binomial likelihood on (ref_reads, alt_reads), so it carries GT, GQ, the
    three PLs and outlier_lod. A row whose reads do not separate the genotypes is
    labelled (FILTER = lowGQ, variant_name = lowqual) rather than deleted, and a
    covered cell that simply does not carry the indel is a called homref rather
    than lowqual. Each called indel allele is genotyped independently (no
    SNV-style per-site collapse), so one row per (barcode, allele) is sufficient.

    hom-ref is per allele:
    n is ref_reads plus this allele's alt_reads, and ref_reads means only
    "carries no indel", since the anchor pileup never records the base. So
    a cell het for a different indel here scores 0/0, and so does one with
    a substitution at the anchor. Other indels and anchor-deleted reads sit
    outside n, below total_reads.

    Depth, VAF, and mean mapping quality all will come along; mean_bq is null
    (indels have no single anchor base quality) so this table unions cleanly with
    the SNV genotype table!

    eps is the SNV substitution error rate applied to an indel, which is
    conservative. indel assay error is the higher of the two, so a real carrier
    cell is never helped by it. It is a flat value here and an indel has no
    anchor base quality to raise it with, so no bq_col is passed.
    """
    return (
        alt_long.with_columns(
            genotype_likelihood_expr(REF_READS, ALT_READS, eps=eps, rho=rho),
            pl.lit(None, dtype=pl.Float64).alias(MEAN_BQ),
        )
        .unnest(_GL)
        .with_columns(genotype_filter_expr(GQ, min_gq).alias(FILTER))
        .with_columns(*genotype_label_exprs(REF_OUT))
        .select(
            BARCODE,
            POS_LOCAL,
            FEATURE,
            REF_OUT,
            ALT_OUT,
            TOTAL_READS,
            REF_READS,
            ALT_READS,
            VAF,
            MEAN_BQ,
            MEAN_MQ,
            GT,
            GQ,
            PL_HOM_REF,
            PL_HET,
            PL_HOM_ALT,
            OUTLIER_LOD,
            ZYGOSITY,
            VARIANT_NAME,
            FILTER,
            INDEL_FEATURE,
        )
        .unique()
    )
