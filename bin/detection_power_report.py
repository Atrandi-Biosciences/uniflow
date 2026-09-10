#!/usr/bin/env python3
"""Per-run detection-power report

Note: Read the detection_power.py script for the source of this code.

Reads the per-cell genotype table this run produced and answers "given what
this run measured, what can this run detect?": callability at the target GQ,
the false-carrier budget, the minimum detectable clone with the allelic-dropout
correction applied, the bulk discovery ceiling in carrier cells, and which one
of those is the binding constraint.

The genotype-likelihood parameters are the ones the run was called with, passed
by conf/variant.config so the report describes the model that was actually
applied. ADO is measured from the data (pseudobulk het sites) and falls back
to the default value with a note when the run cannot measure it.

SNV only. The indel genotype table borrows the SNV eps , so a separate indel
power number would restate this one with an error rate known to be wrong for indels.

Outputs, both to dev/:
  * <sample>_detection_power.txt   the report
  * <sample>_detection_power.csv   the same numbers, one row per quantity
"""

import argparse

import polars as pl
from lib.common_const import BARCODE
from lib.modality.dna.detection_power import (
    DEFAULT_BULK_MIN_ALT_READS,
    DEFAULT_BULK_MIN_VAF,
    DEFAULT_FAMILY_ALPHA,
    DEFAULT_GQ_CALIBRATION,
    DEFAULT_MIN_CALLERS,
    GQ_CALIBRATION_CHOICES,
    compute_power,
    render_text,
    to_table,
)
from lib.modality.dna.processing import (
    ALT_READS,
    DEFAULT_ERROR_EPS,
    DEFAULT_ERROR_RHO,
    DEFAULT_MIN_GQ,
    FEATURE,
    GQ,
    N_CALLERS,
    POS_LOCAL,
    REF_READS,
)

parser = argparse.ArgumentParser(
    description="Per-run detection-power report (ADR-F14)."
)
parser.add_argument("variants_count_parquet")
parser.add_argument("sample_name")
parser.add_argument("--min-gq", type=float, default=DEFAULT_MIN_GQ)
parser.add_argument("--error-eps", type=float, default=DEFAULT_ERROR_EPS)
parser.add_argument("--error-rho", type=float, default=DEFAULT_ERROR_RHO)
parser.add_argument(
    "--bulk-min-vaf",
    type=float,
    default=DEFAULT_BULK_MIN_VAF,
    help="Bulk discovery VAF floor, to express as carrier cells.",
)
parser.add_argument(
    "--bulk-min-alt-reads",
    type=int,
    default=DEFAULT_BULK_MIN_ALT_READS,
    help="Bulk discovery pooled AD[ALT] floor, to express as carrier cells.",
)
parser.add_argument(
    "--family-alpha",
    type=float,
    default=DEFAULT_FAMILY_ALPHA,
    help="Expected number of sites panel-wide where noise alone fakes a clone.",
)
parser.add_argument(
    "--min-callers",
    type=int,
    default=DEFAULT_MIN_CALLERS,
    help="Caller support defining the second reporting scope for the budget.",
)
parser.add_argument(
    "--ado",
    type=float,
    default=None,
    help="Override the measured ADO rate instead of estimating it from this run.",
)
parser.add_argument(
    "--gq-calibration",
    choices=GQ_CALIBRATION_CHOICES,
    default=DEFAULT_GQ_CALIBRATION,
    help=(
        "How rho, and therefore GQ, was calibrated for this sample. rho cannot be "
        "estimated without a truth set, so a sample that has none inherits the "
        "existing constant and its GQ is optimistic by however much this sample is "
        "more overdispersed than the calibration run. Note thtat this not a "
        "measurement, thus the data cannot settle it."
    ),
)
args = parser.parse_args()

scan = pl.scan_parquet(args.variants_count_parquet)
available = scan.collect_schema().names()

# ref_reads is the column the whole report is built on: n = ref + alt,
# never total_reads. A table without it predates the genotype
# likelihood and cannot be given one retroactively, since the third-allele
# counts it excludes are not recoverable from the row.
missing = [c for c in (REF_READS, ALT_READS) if c not in available]
if missing:
    raise SystemExit(
        f"{args.variants_count_parquet} has no {', '.join(missing)} column: it predates "
        "the per-cell genotype likelihood, so its depth cannot be reconstructed."
    )

# Only the columns the report reads, so a production-scale genotype table only takes
# a fraction of its width in memory.
# GQ is optional: it is reported as the observed callable fraction and never fed back
# into a derived number.
wanted = [BARCODE, FEATURE, POS_LOCAL, REF_READS, ALT_READS]
optional = [c for c in (GQ, N_CALLERS) if c in available]
frame = scan.select(wanted + optional).collect()

report = compute_power(
    frame,
    sample=args.sample_name,
    eps=args.error_eps,
    rho=args.error_rho,
    min_gq=args.min_gq,
    bulk_min_vaf=args.bulk_min_vaf,
    bulk_min_alt_reads=args.bulk_min_alt_reads,
    family_alpha=args.family_alpha,
    min_callers=args.min_callers,
    ado_override=args.ado,
    gq_calibration=args.gq_calibration,
)

text = render_text(report)
print(text)
with open(f"{args.sample_name}_detection_power.txt", "w") as handle:
    handle.write(text)
to_table(report).write_csv(f"{args.sample_name}_detection_power.csv")
