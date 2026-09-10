# lib/modality/dna/detection_power.py
"""Detection power is basically what can this run detect?

Note: This is heavily inspired from the previous work on the target enrichment sequencing
of whole genome. While it is not the same technology, but the nature of data, variant calling,
zygosity, and clones exist as a concept. The only new additions are allelic dropout concept
which are amplicon focused only. This report is highly experimental and require extensive validation,
verification, and fine tuning. All the numbers here are derived from investigating a 10M downsampled
experiment.

This module asnwer the question:
"Given the depth this run actually has, the cells it actually called and the
dropout it actually shows, what is the smallest clone it could have seen,
and what is stopping it from seeing a smaller one?"

Everything here reads the genotype table the run just produced
(variants_count.parquet) and the genotype-likelihood parameters that produced
it, so the report describes the model that was actually applied rather than an
assumed one.

Of course nothing is hardcoded per assay. The two quantities that cannot be
derived from an untruthed run (eps, rho) arrive as arguments, and ADO is
measured from the data by the pseudobulk discovery phase route with a fallback.

The six outputs are, in the order they require and drive each other:

  1. callability: can the provided depth reach the target GQ at all?
  2. false carriers: eps x cells x sites, the limits that set min_gq
  3. C_min: Poisson counting against that limit
  4. the ADO correction: dropout removes carrier cells for that variant, it does not add noise
  5. the bulk discovery phase ceiling: bulk_min_vaf / AD floor expressed in carrier cells
  6. the binding constraint, named

Read section 5 of docs/adr/variant_filtering.md before changing any of it. The
appendix (section 7) defines every term used here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

import numpy as np
import polars as pl
from scipy.special import gammaln

from lib.common_const import BARCODE
from lib.modality.dna.processing import (
    ALT_READS,
    FEATURE,
    GQ,
    N_CALLERS,
    POS_LOCAL,
    REF_READS,
    THETA_HET,
    _beta_binomial_log_pmf,
    genotype_likelihood_arrays,
)

# Pseudobulk germline-het band. A site whose pooled VAF over
# all cells end up here is a germline heterozygote in this sample; assay's seq noise
# level sits at ~1e-3 and a hom-alt at ~1, so the range is wide but unambiguous.
# It is not a genotype call, only a site selector for the ADO estimator.
HET_SITE_BAND: Final[tuple[float, float]] = (0.35, 0.65)

# A dropout event: a cell at a germline het site whose VAF collapsed to one
# allele. Although this number is meaningless without proper testing
# (0.05 / 0.10 / 0.20 give ADO 9.5% / 14.0% / 21.5% on a downsampled 10M run).
ADO_BAND: Final[float] = 0.10

# A cell whose depth cannot resolve an allele fraction cannot report a dropout
# either, so the estimator only looks at cells with this much ref+alt evidence.
ADO_MIN_DP: Final[int] = 20

# Below these the pseudobulk estimate is noise, and the configured fallback is
# used instead (declared as such in the report; ADO is never just assumed).
ADO_MIN_SITES: Final[int] = 5
ADO_MIN_OBS: Final[int] = 100

# A measured ADO on the 10M run returned 0.14 (see above; the numbers for ADO_BAND),
# the fallback when this run cannot measure its own. Like eps and rho it is a seed,
# not a constant.
DEFAULT_ADO: Final[float] = 0.140

# Family-wise limit for the carrier cell-count test: it is the expected number of sites
# panel-wide at which noise alone would produce a "clone". Split over the sites
# actually tested (Bonferroni), which is what makes C_min scale with the panel.
DEFAULT_FAMILY_ALPHA: Final[float] = 0.05

# The second reporting scope: sites backed by at least this many callers. 2 is
# the effective authenticity filter measured in the false-positive analysis
# the bulk floor is deliberately permissive, so the all-sites budget and this
# one answer different questions.
DEFAULT_MIN_CALLERS: Final[int] = 2

# The bulk discovery floor, mirrored from conf/variant.config so the report can
# express it in carrier cells. Overridden from the CLI by the pipeline.
DEFAULT_BULK_MIN_VAF: Final[float] = 1e-4
DEFAULT_BULK_MIN_ALT_READS: Final[int] = 20

# A run below this callable fraction is depth-limited, and no amount of
# threshold tuning recovers it.
CALLABLE_TARGET: Final[float] = 0.90

DEPTH: Final[str] = "n"
GT_HOM_REF_INDEX: Final[int] = 0

# Thus rho is a constant calibrated once on a truth VCF run, and this says
# whether the run reading it is that run. It is not something to measure and thus
# the data cannot settle on its own, so it is the one value that is a fact about the
# SAMPLE rather than a tuning choice.
# That is why it arrives as the gq_calibration pipeline parameter rather than as a
# conf/variant.config threshold like every other value here.
GQ_TRUTH_CALIBRATED: Final[str] = "truth-calibrated"
GQ_INHERITED: Final[str] = "inherited"
GQ_CALIBRATION_CHOICES: Final[tuple[str, ...]] = (
    GQ_TRUTH_CALIBRATED,
    GQ_INHERITED,
)
DEFAULT_GQ_CALIBRATION: Final[str] = GQ_INHERITED

CONSTRAINT_BULK: Final[str] = "bulk discovery floor"
CONSTRAINT_EPS: Final[str] = "eps (false carrier-cell noise at the chosen GQ)"
CONSTRAINT_DEPTH: Final[str] = "per-cell depth"
CONSTRAINT_ADO: Final[str] = "allelic dropout"
CONSTRAINT_POISSON: Final[str] = "Poisson counting"


def log_choose(n, k):
    """log C(n, k). The likelihood drops this term because it cancels between
    diploid states; a probability over k does not get to drop it."""
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)


def state_pmf(n: int, theta: float, rho: float) -> np.ndarray:
    """P(k ALT of n) for k = 0..n under one genotype state. Sums to 1.

    theta 0 or 1 is a point mass. The beta-binomial's alpha goes to 0 there
    and betaln(0, .) is infinite, so the general expression returns NaN for
    the one outcome that has all the mass. An error-free model (eps = 0) is
    also possible, and quite valid thing to ask this for!
    """
    pmf = np.zeros(n + 1, dtype=np.float64)
    if theta <= 0.0:
        pmf[0] = 1.0
        return pmf
    if theta >= 1.0:
        pmf[n] = 1.0
        return pmf
    k = np.arange(n + 1, dtype=np.float64)
    return np.exp(
        log_choose(float(n), k) + _beta_binomial_log_pmf(k, float(n), theta, rho)
    )


def _calls_at_depth(n: int, eps: float, rho: float, theta_het: float):
    """(gq, gt_index) for every possible k at depth n; the call the pipeline
    would make for each outcome, from the same function the pipeline uses."""
    k = np.arange(n + 1, dtype=np.int64)
    _, _, _, gq, gt_index, _ = genotype_likelihood_arrays(
        n - k, k, eps=eps, rho=rho, theta_het=theta_het
    )
    return gq, gt_index


def max_gq_by_state(
    n: int, eps: float, rho: float, theta_het: float = THETA_HET
) -> tuple[float, float, float]:
    """The best GQ each genotype state can reach at depth n (0 if unreachable).

    The cieling for this is 1 and it is computed rather than approximated
    approximated by the '~3 x DP' rule: at rho = 0.091 a het saturates near
    51 whatever the depth, while a hom-ref keeps climbing, so one rule cannot
    describe both arms.
    """
    if n <= 0:
        return (0.0, 0.0, 0.0)
    gq, gt_index = _calls_at_depth(n, eps, rho, theta_het)
    return tuple(
        float(gq[gt_index == state].max()) if (gt_index == state).any() else 0.0
        for state in (0, 1, 2)
    )


def depth_for_gq(
    target_gq: float,
    eps: float,
    rho: float,
    theta_het: float = THETA_HET,
    max_depth: int = 5000,
) -> tuple[int | None, int | None]:
    """Smallest depth at which a clean hom-ref / a perfect het can reach
    target_gq. None when the state saturates below the target at any depth,
    which is the het arm's situation at any realistic rho."""
    reached: list[int | None] = [None, None]
    n = 1
    while n <= max_depth and (reached[0] is None or reached[1] is None):
        ceiling = max_gq_by_state(n, eps, rho, theta_het)
        for i, state in enumerate((0, 1)):
            if reached[i] is None and ceiling[state] >= target_gq:
                reached[i] = n
        n = n * 2 if n >= 64 else n + 1
    # The doubling above overshoots so we walk back to the exact smallest depth.
    for i, state in enumerate((0, 1)):
        if reached[i] is not None:
            while (
                reached[i] > 1
                and max_gq_by_state(reached[i] - 1, eps, rho, theta_het)[state]
                >= target_gq
            ):
                reached[i] -= 1
    return reached[0], reached[1]


def false_carrier_rate(
    n: int, eps: float, rho: float, min_gq: float, theta_het: float = THETA_HET
) -> float:
    """P(a TRUE hom-ref cell at depth n is called non-reference at GQ >= min_gq).

    Exact: enumerate every k, ask the pipeline's own likelihood what it would
    call, and weight by the probability of that k under the error model. No
    simulation and no normal approximation, both of which fail in exactly the
    tail this number lives in.
    """
    if n <= 0:
        return 0.0
    gq, gt_index = _calls_at_depth(n, eps, rho, theta_het)
    pmf = state_pmf(n, eps, rho)
    return float(pmf[(gt_index != GT_HOM_REF_INDEX) & (gq >= min_gq)].sum())


def poisson_upper_tail(c: int, lam: float) -> float:
    """P(X >= c) for X ~ Poisson(lam), summed exactly in log space."""
    if c <= 0:
        return 1.0
    if lam <= 0.0:
        return 0.0
    i = np.arange(c, dtype=np.float64)
    below = np.exp(-lam + i * np.log(lam) - gammaln(i + 1.0)).sum()
    return float(max(0.0, 1.0 - below))


def poisson_min_count(lam: float, alpha: float, max_count: int = 10_000) -> int:
    """Smallest carrier cell count whose upper tail under noise alone is <= alpha.

    Never below 1: a single carrier cell in a run with no false-carrier
    expectation is already significant, and reporting C_min = 0 would say the
    assay detects clones it never observed.
    """
    c = 1
    while c < max_count and poisson_upper_tail(c, lam) > alpha:
        c += 1
    return c


def measure_ado(
    frame: pl.DataFrame,
    band: float = ADO_BAND,
    min_dp: int = ADO_MIN_DP,
    het_site_band: tuple[float, float] = HET_SITE_BAND,
) -> tuple[float | None, int, int]:
    """ADO from pseudobulk germline het sites; (rate, n_sites, n_observations).

    A site whose pooled VAF over all cells is ~0.5 is a germline het, so every adequately
    covered cell at that site should show both alleles, and the ones that do not
    dropped one. Weaker than a supplied donor VCF (it is confounded in a donor
    mixture, where a site het in one donor and hom-ref in the other also pools
    near 0.25) and it is deliberately not used to genotype anything.

    Returns rate None when the run has too few het sites or too little depth to
    say; the caller then declares the fallback rather than reporting a number
    built on three cells.
    """
    if frame.height == 0:
        return (None, 0, 0)

    per_site = frame.group_by(FEATURE, POS_LOCAL).agg(
        pl.col(REF_READS).sum().alias("site_ref"),
        pl.col(ALT_READS).sum().alias("site_alt"),
    )
    het_sites = per_site.filter((pl.col("site_ref") + pl.col("site_alt")) > 0).filter(
        (pl.col("site_alt") / (pl.col("site_ref") + pl.col("site_alt"))).is_between(
            *het_site_band
        )
    )
    if het_sites.height == 0:
        return (None, 0, 0)

    observations = (
        frame.join(het_sites.select(FEATURE, POS_LOCAL), on=[FEATURE, POS_LOCAL])
        .with_columns((pl.col(REF_READS) + pl.col(ALT_READS)).alias(DEPTH))
        .filter(pl.col(DEPTH) >= min_dp)
    )
    if observations.height < ADO_MIN_OBS or het_sites.height < ADO_MIN_SITES:
        return (None, het_sites.height, observations.height)

    cell_vaf = (observations[ALT_READS] / observations[DEPTH]).to_numpy()
    rate = float(((cell_vaf <= band) | (cell_vaf >= 1.0 - band)).mean())
    return (rate, het_sites.height, observations.height)


@dataclass
class PowerReport:
    """Every number the report prints, so the text and the CSV cannot diverge."""

    sample: str
    n_cells: int
    n_sites: int
    n_observations: int
    # parameters as applied
    eps: float
    rho: float
    gq_calibration: str
    theta_het: float
    min_gq: float
    # 1. callability
    depth_quantiles: dict[str, float]
    dp_for_gq_homref: int | None
    dp_for_gq_het: int | None
    callable_fraction_homref: float
    callable_fraction_het: float
    observed_pass_fraction: float | None
    het_gq_ceiling: float
    # 2. false carrier cells
    false_carrier_rate_mean: float
    expected_false_carriers_panel: float
    lambda_site_median: float
    lambda_site_max: float
    # 3. counting, over the whole panel and over the authenticity-filtered ones
    family_alpha: float
    alpha_per_site: float
    c_min: int
    min_callers: int
    scoped_n_sites: int | None
    scoped_alpha_per_site: float | None
    scoped_c_min: int | None
    scoped_expected_false_carriers: float | None
    # 4. ADO Allelic dropout
    ado: float
    ado_measured: bool
    ado_sites: int
    ado_observations: int
    c_true_counting: float
    clonal_fraction_counting: float
    # 5. bulk ceiling
    bulk_min_vaf: float
    bulk_min_alt_reads: int
    c_bulk_vaf: float
    c_bulk_alt_reads: float
    c_true_bulk: float
    # 6. the answer
    c_true_min: float
    clonal_fraction_min: float
    binding_constraint: str
    notes: list[str] = field(default_factory=list)


def compute_power(
    frame: pl.DataFrame,
    sample: str,
    eps: float,
    rho: float,
    min_gq: float,
    theta_het: float = THETA_HET,
    bulk_min_vaf: float = DEFAULT_BULK_MIN_VAF,
    bulk_min_alt_reads: int = DEFAULT_BULK_MIN_ALT_READS,
    family_alpha: float = DEFAULT_FAMILY_ALPHA,
    min_callers: int = DEFAULT_MIN_CALLERS,
    ado_override: float | None = None,
    gq_calibration: str = DEFAULT_GQ_CALIBRATION,
) -> PowerReport:
    """The whole report, from one genotype table plus the parameters that built it."""
    notes: list[str] = []
    if gq_calibration not in GQ_CALIBRATION_CHOICES:
        raise ValueError(
            f"gq_calibration must be one of {GQ_CALIBRATION_CHOICES}, "
            f"got {gq_calibration!r}. It is a declaration about how rho was "
            "calibrated and cannot be guessed from the data."
        )
    if gq_calibration == GQ_INHERITED:
        notes.append(
            "GQ is 'inherited' from an existing reference run, not calibrated on this one: "
            "rho cannot be estimated without a truth set. GT and outlier_lod are unaffected"
            ", but every GQ here is optimistics if this sample is more overdispersed than "
            "the calibration run, and so are the het ceiling, the callable fractions and the"
            "depth targets derived from it."
        )
    work = frame.with_columns(
        (pl.col(REF_READS).fill_null(0) + pl.col(ALT_READS).fill_null(0)).alias(DEPTH)
    )
    n_cells = work[BARCODE].n_unique()
    n_sites = work.select(FEATURE, POS_LOCAL).n_unique()
    n_obs = work.height

    depths = work[DEPTH].to_numpy()
    quantiles = {
        f"p{int(q * 100)}": float(np.quantile(depths, q))
        for q in (0.1, 0.25, 0.5, 0.75, 0.9)
    }
    quantiles["mean"] = float(depths.mean())
    quantiles["max"] = float(depths.max())

    # 1. callability. Per distinct depth, since the ceiling is a function of
    # depth alone then a handful of distinct depths stands behind every observation.
    distinct, inverse, counts = np.unique(
        depths, return_inverse=True, return_counts=True
    )
    ceilings = np.array(
        [max_gq_by_state(int(n), eps, rho, theta_het) for n in distinct],
        dtype=np.float64,
    )
    weight = counts / counts.sum()
    callable_homref = float(weight[ceilings[:, 0] >= min_gq].sum())
    callable_het = float(weight[ceilings[:, 1] >= min_gq].sum())
    dp_homref, dp_het = depth_for_gq(min_gq, eps, rho, theta_het)
    het_ceiling = float(max_gq_by_state(int(distinct.max()), eps, rho, theta_het)[1])

    observed_pass = None
    if GQ in work.columns:
        observed_pass = float((work[GQ].fill_null(0.0).to_numpy() >= min_gq).mean())

    if callable_het == 0.0:
        notes.append(
            "NO depth in this run can reach the target GQ as a heterozygote; "
            "every carrier call is lowqual by arithmetic, not by evidence."
        )

    # 2. false carrier cells. Rate per distinct depth, weighted by how often that
    # depth occurs, then multiplied out over the observations that exist.
    rates = np.array(
        [false_carrier_rate(int(n), eps, rho, min_gq, theta_het) for n in distinct]
    )
    rate_per_obs = rates[inverse]
    false_rate_mean = float(rate_per_obs.mean())
    expected_false = float(rate_per_obs.sum())

    per_site_rate = work.with_columns(pl.Series("_rate", rate_per_obs))
    site_lambda = (
        per_site_rate.group_by(FEATURE, POS_LOCAL)
        .agg(pl.col("_rate").sum().alias("lam"))["lam"]
        .to_numpy()
    )
    lam_median = float(np.median(site_lambda)) if site_lambda.size else 0.0
    lam_max = float(site_lambda.max()) if site_lambda.size else 0.0

    # 3. counting, at the worst site rather than the median one: the panel is
    # reported as a whole, so the site most likely to show a clone sets
    # the threshold the whole panel is read at.
    alpha_per_site = family_alpha / max(n_sites, 1)
    c_min = poisson_min_count(lam_max, alpha_per_site)

    # The same limit over the authenticity-filtered ones, when the table
    # carries caller support. Measured on the 10M run: cutting 2,071 sites to
    # 130 loosens alpha 16x and does not move C_min (the Poisson tail is too
    # steep for that), while it cuts expected false carrier cells 48 -> 3. So the
    # scope is a statement about how much noise a reader sees, not about
    # sensitivity, and printing both is what makes that visible. A site counts
    # as supported if ANY of its rows is, because a hom-ref cell at a called
    # site carries a null ALT and so no caller support of its own.
    scoped_sites = scoped_alpha = scoped_c_min = scoped_false = None
    if N_CALLERS in work.columns:
        supported = (
            per_site_rate.group_by(FEATURE, POS_LOCAL)
            .agg(
                pl.col("_rate").sum().alias("lam"),
                pl.col(N_CALLERS).fill_null(0).max().alias("support"),
            )
            .filter(pl.col("support") >= min_callers)
        )
        if supported.height > 0:
            scoped_lam = supported["lam"].to_numpy()
            scoped_sites = supported.height
            scoped_alpha = family_alpha / scoped_sites
            scoped_c_min = poisson_min_count(float(scoped_lam.max()), scoped_alpha)
            scoped_false = float(scoped_lam.sum())

    # 4. ADO. Dropout removes carrier cells, so it inflates every detection limit.
    if ado_override is not None:
        ado, ado_measured, ado_sites, ado_obs = (ado_override, False, 0, 0)
        notes.append(
            f"ADO was supplied as {ado_override:.3f}, not measured on this run."
        )
    else:
        measured, ado_sites, ado_obs = measure_ado(work)
        if measured is None:
            ado, ado_measured = DEFAULT_ADO, False
            notes.append(
                f"ADO not measurable on this run ({ado_sites} pseudobulk het sites, "
                f"{ado_obs} cell observations at DP >= {ADO_MIN_DP}); falling back to "
                f"the default value {DEFAULT_ADO:.3f}. Every clone size below "
                "inherits that assumption."
            )
        else:
            ado, ado_measured = measured, True

    # A true clone of C cells presents as C x (1 - ADO) carriers, of which only
    # the callable fraction can be genotyped at the target GQ. Both shrink the
    # observed count, so both inflate the true clone needed to produce C_min.
    observable = (1.0 - ado) * callable_het
    c_true_counting = float(c_min / observable) if observable > 0 else float("inf")

    # 5. the bulk ceiling. A variant the bulk layer never discovers is never
    # piled up per cell, so this is a hard ceiling on the per-cell layer whatever
    # its own power. Two floors apply; the binding one is the larger.
    mean_depth = float(depths.mean()) if depths.size else 0.0
    c_bulk_vaf = bulk_min_vaf * n_cells / theta_het
    c_bulk_ad = (
        bulk_min_alt_reads / (theta_het * mean_depth)
        if mean_depth > 0
        else float("inf")
    )
    # Total dropout leaves no ALT read anywhere, in bulk or per cell, so the
    # limit is infinite rather than a division by zero.
    c_true_bulk = (
        float(max(c_bulk_vaf, c_bulk_ad) / (1.0 - ado)) if ado < 1.0 else float("inf")
    )

    # 6. the binding constraint, by an ordered we can compare two runs
    c_true_min = max(c_true_counting, c_true_bulk)
    if c_true_bulk >= c_true_counting:
        binding = CONSTRAINT_BULK
    elif c_min > 1:
        binding = CONSTRAINT_EPS
    elif callable_het < CALLABLE_TARGET:
        binding = CONSTRAINT_DEPTH
    elif ado > ADO_BAND:
        binding = CONSTRAINT_ADO
    else:
        binding = CONSTRAINT_POISSON

    return PowerReport(
        sample=sample,
        n_cells=n_cells,
        n_sites=n_sites,
        n_observations=n_obs,
        eps=eps,
        rho=rho,
        gq_calibration=gq_calibration,
        theta_het=theta_het,
        min_gq=min_gq,
        depth_quantiles=quantiles,
        dp_for_gq_homref=dp_homref,
        dp_for_gq_het=dp_het,
        callable_fraction_homref=callable_homref,
        callable_fraction_het=callable_het,
        observed_pass_fraction=observed_pass,
        het_gq_ceiling=het_ceiling,
        false_carrier_rate_mean=false_rate_mean,
        expected_false_carriers_panel=expected_false,
        lambda_site_median=lam_median,
        lambda_site_max=lam_max,
        family_alpha=family_alpha,
        alpha_per_site=alpha_per_site,
        c_min=c_min,
        min_callers=min_callers,
        scoped_n_sites=scoped_sites,
        scoped_alpha_per_site=scoped_alpha,
        scoped_c_min=scoped_c_min,
        scoped_expected_false_carriers=scoped_false,
        ado=ado,
        ado_measured=ado_measured,
        ado_sites=ado_sites,
        ado_observations=ado_obs,
        c_true_counting=c_true_counting,
        clonal_fraction_counting=c_true_counting / n_cells if n_cells else float("inf"),
        bulk_min_vaf=bulk_min_vaf,
        bulk_min_alt_reads=bulk_min_alt_reads,
        c_bulk_vaf=c_bulk_vaf,
        c_bulk_alt_reads=c_bulk_ad,
        c_true_bulk=c_true_bulk,
        c_true_min=c_true_min,
        clonal_fraction_min=c_true_min / n_cells if n_cells else float("inf"),
        binding_constraint=binding,
        notes=notes,
    )


def _fmt(value: float | None, spec: str = ".4g") -> str:
    if value is None:
        return "unreachable"
    if isinstance(value, float) and not np.isfinite(value):
        return "inf"
    return format(value, spec)


# The report is one template so its layout is literal. It can move to a
# html report later on. For now, it is just a template to be filled!
_REPORT = """\
Detection power for {r.sample}
{rule}

  cells {r.n_cells}   sites {r.n_sites}   cell-site observations {r.n_observations}
  genotype model as applied: eps {r.eps:.3g}, rho {r.rho:.3g}, theta_het {r.theta_het:.3g}, min_gq {r.min_gq:.0f}
  GQ calibration: {r.gq_calibration} -- {gq_provenance}

1. CALLABILITY  (can the observed depth reach the target GQ at all)
   depth n = ref+alt per cell-site:  {depths}
   depth needed for GQ {r.min_gq:.0f}:  hom-ref {dp_homref}   het {dp_het}
   observations that can reach it:   hom-ref {r.callable_fraction_homref:>7.2%}   het {r.callable_fraction_het:>7.2%}
   het GQ ceiling at this run's deepest cell-site: {r.het_gq_ceiling:.1f} (rho caps the het arm; the hom-ref arm does not saturate){observed}

2. FALSE-CARRIER LIMIT  (a true hom-ref cell called non-reference)
   per-cell-site rate, depth-weighted: {r.false_carrier_rate_mean:.3g}
   expected false carriers panel-wide: {r.expected_false_carriers_panel:.2f}
   per site: median {r.lambda_site_median:.3g}, worst {r.lambda_site_max:.3g}

3. MINIMUM DETECTABLE CARRIER CELL COUNT
   family-wise alpha {r.family_alpha:.3g} over {r.n_sites} sites -> alpha per site {r.alpha_per_site:.3g}
   C_min (observed carrier cells at the worst site): {r.c_min}{scoped}

4. THE ADO CORRECTION  (dropout removes carrier cells, it does not add noise)
   ADO {r.ado:.3f} {ado_source}
   observable fraction of a true clone: (1 - ADO) x callable_het = {observable:.3f}
   minimum detectable true clone, counting arm: {c_counting} cells = {pct_counting}% of cells

5. THE BULK DISCOVERY CEILING  (a variant never discovered is never counted)
   bulk_min_vaf {r.bulk_min_vaf:.3g} at {r.n_cells} cells -> {c_bulk_vaf} carrier cells
   AD[ALT] >= {r.bulk_min_alt_reads} at mean depth {mean_depth} -> {c_bulk_ad} carrier cells
   minimum detectable TRUE clone, bulk arm: {c_bulk} cells

6. THE BINDING CONSTRAINT
   >>> {r.binding_constraint} <<<
   minimum detectable clone for this run: {c_min_true} cells = {pct_min}% clonal fraction{notes}
"""

_OBSERVED = """
   observations that did reach it (gq >= min_gq): {r.observed_pass_fraction:.2%}"""

_SCOPED = """
   same limitation over sites with >= {r.min_callers} callers: {r.scoped_n_sites} sites -> alpha per site {r.scoped_alpha_per_site:.3g}, C_min {r.scoped_c_min}, expected false carrier cells {r.scoped_expected_false_carriers:.2f}
   (site mainly moves the false carrier cell burden : the Poisson tail
    is steep, so a 16x looser alpha typically leaves C_min where it is.
    Filtering sites buys cleanliness, not sensitivity.)"""

_ADO_MEASURED = "MEASURED on {r.ado_sites} pseudobulk het sites, {r.ado_observations} cell observations at DP >= {min_dp}, band {band:.2f}"

# GQ is only accurate on a run that has a truth set, so the report says
# which case it is in rather than leaving the reader to assume the good one.
_GQ_PROVENANCE = {
    GQ_TRUTH_CALIBRATED: (
        "rho was calibrated on this run's own truth set, so GQ is accurate here "
        "and this run is a calibration source for others"
    ),
    GQ_INHERITED: (
        "rho comes from a reference run, so GQ is optimistic if this sample is "
        "more overdispersed than that one. GT and outlier_lod are unaffected"
    ),
}

# (labels in the report with depth_quantiles)
_DEPTH_FIELDS = (
    ("p10", "p10"),
    ("p25", "p25"),
    ("median", "p50"),
    ("p75", "p75"),
    ("p90", "p90"),
    ("max", "max"),
)


def render_text(report: PowerReport) -> str:
    """The human-readable report. Every number is also in the CSV."""
    r = report
    q = r.depth_quantiles
    ado_source = (
        _ADO_MEASURED.format(r=r, min_dp=ADO_MIN_DP, band=ADO_BAND)
        if r.ado_measured
        else "ASSUMED (see notes)"
    )
    return _REPORT.format(
        r=r,
        rule="=" * 72,
        depths="  ".join(f"{label} {_fmt(q[key])}" for label, key in _DEPTH_FIELDS),
        dp_homref=_fmt(r.dp_for_gq_homref),
        dp_het=_fmt(r.dp_for_gq_het),
        observed=_OBSERVED.format(r=r) if r.observed_pass_fraction is not None else "",
        scoped=_SCOPED.format(r=r) if r.scoped_n_sites is not None else "",
        ado_source=ado_source,
        gq_provenance=_GQ_PROVENANCE[r.gq_calibration],
        observable=(1.0 - r.ado) * r.callable_fraction_het,
        c_counting=_fmt(r.c_true_counting, ".1f"),
        pct_counting=_fmt(r.clonal_fraction_counting * 100, ".3g"),
        mean_depth=_fmt(q["mean"]),
        c_bulk_vaf=_fmt(r.c_bulk_vaf, ".1f"),
        c_bulk_ad=_fmt(r.c_bulk_alt_reads, ".1f"),
        c_bulk=_fmt(r.c_true_bulk, ".1f"),
        c_min_true=_fmt(r.c_true_min, ".1f"),
        pct_min=_fmt(r.clonal_fraction_min * 100, ".3g"),
        notes=(
            "\n\nNOTES\n" + "\n".join(f"   - {n}" for n in r.notes) if r.notes else ""
        ),
    )


def to_table(report: PowerReport) -> pl.DataFrame:
    """One row per quantity: sample, quantity, value, unit, provenance.

    Long rather than wide so a cross-run store can append runs without a schema
    migration every time a quantity is added (see docs/results_accumulation_design.md).
    """
    r = report
    ado_provenance = "measured" if r.ado_measured else "assumed"
    rows: list[tuple[str, float | None, str, str]] = [
        ("n_cells", r.n_cells, "cells", "measured"),
        ("n_sites", r.n_sites, "sites", "measured"),
        ("n_observations", r.n_observations, "cell-sites", "measured"),
        ("eps", r.eps, "fraction", "configured"),
        # rho has its calibration source rather than a flat
        # "configured", since that is what tells a cross-run store whether this
        # run's GQ is accurate or inherited.
        ("rho", r.rho, "fraction", r.gq_calibration),
        ("theta_het", r.theta_het, "fraction", "fixed"),
        ("min_gq", r.min_gq, "phred", "configured"),
        ("depth_median", r.depth_quantiles["p50"], "reads", "measured"),
        ("depth_mean", r.depth_quantiles["mean"], "reads", "measured"),
        ("depth_p10", r.depth_quantiles["p10"], "reads", "measured"),
        ("depth_p90", r.depth_quantiles["p90"], "reads", "measured"),
        ("depth_for_target_gq_homref", r.dp_for_gq_homref, "reads", "derived"),
        ("depth_for_target_gq_het", r.dp_for_gq_het, "reads", "derived"),
        ("callable_fraction_homref", r.callable_fraction_homref, "fraction", "derived"),
        ("callable_fraction_het", r.callable_fraction_het, "fraction", "derived"),
        ("observed_pass_fraction", r.observed_pass_fraction, "fraction", "measured"),
        ("het_gq_ceiling", r.het_gq_ceiling, "phred", "derived"),
        ("false_carrier_rate", r.false_carrier_rate_mean, "per cell-site", "derived"),
        (
            "expected_false_carriers",
            r.expected_false_carriers_panel,
            "cells",
            "derived",
        ),
        ("false_carriers_per_site_median", r.lambda_site_median, "cells", "derived"),
        ("false_carriers_per_site_max", r.lambda_site_max, "cells", "derived"),
        ("family_alpha", r.family_alpha, "probability", "configured"),
        ("alpha_per_site", r.alpha_per_site, "probability", "derived"),
        ("c_min_observed", r.c_min, "cells", "derived"),
        ("min_callers", r.min_callers, "callers", "configured"),
        ("scoped_n_sites", r.scoped_n_sites, "sites", "measured"),
        ("scoped_alpha_per_site", r.scoped_alpha_per_site, "probability", "derived"),
        ("scoped_c_min_observed", r.scoped_c_min, "cells", "derived"),
        (
            "scoped_expected_false_carriers",
            r.scoped_expected_false_carriers,
            "cells",
            "derived",
        ),
        ("ado", r.ado, "fraction", ado_provenance),
        ("ado_sites", r.ado_sites, "sites", "measured"),
        ("ado_observations", r.ado_observations, "cell-sites", "measured"),
        ("c_true_counting", r.c_true_counting, "cells", "derived"),
        ("clonal_fraction_counting", r.clonal_fraction_counting, "fraction", "derived"),
        ("bulk_min_vaf", r.bulk_min_vaf, "fraction", "configured"),
        ("bulk_min_alt_reads", r.bulk_min_alt_reads, "reads", "configured"),
        ("c_bulk_vaf", r.c_bulk_vaf, "cells", "derived"),
        ("c_bulk_alt_reads", r.c_bulk_alt_reads, "cells", "derived"),
        ("c_true_bulk", r.c_true_bulk, "cells", "derived"),
        ("c_true_min", r.c_true_min, "cells", "derived"),
        ("clonal_fraction_min", r.clonal_fraction_min, "fraction", "derived"),
    ]
    return pl.DataFrame(
        {
            "sample": [r.sample] * (len(rows) + 1),
            "quantity": [row[0] for row in rows] + ["binding_constraint"],
            "value": [None if row[1] is None else float(row[1]) for row in rows]
            + [None],
            "text": [None] * len(rows) + [r.binding_constraint],
            "unit": [row[2] for row in rows] + ["name"],
            "provenance": [row[3] for row in rows] + ["derived"],
        },
        schema={
            "sample": pl.String,
            "quantity": pl.String,
            "value": pl.Float64,
            "text": pl.String,
            "unit": pl.String,
            "provenance": pl.String,
        },
    )
