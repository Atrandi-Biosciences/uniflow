#!/usr/bin/env python3
"""Convert a merged consensus VCF into a long-format parquet table.

One row per (sample, caller, variant) where the caller has a non-missing
FORMAT/DP at that record. All per-caller signal is read out of the merged
VCF's FORMAT columns (AD, DP, VAF, VARQUAL, FILTER); caller-specific INFO
fields are left empty in v1.

Usage:
    vcf_to_parquet.py <merged.vcf.gz> <sample_name> <output.parquet>
"""

import sys

import polars as pl
import pysam

SCHEMA = {
    "sample_name": pl.String,
    "caller": pl.String,
    "chrom": pl.String,
    "pos": pl.Int64,
    "ref": pl.String,
    "alt": pl.String,
    "type": pl.String,
    "qual": pl.Float64,
    "dp": pl.Int64,
    "ad_ref": pl.Int64,
    "ad_alt": pl.Int64,
    "af": pl.Float64,
    "caller_specific": pl.String,
}


def _variant_type(ref, alt):
    """Classify a normalized (split + left-aligned) REF/ALT pair.

    BCFTOOLS_NORM splits multi-allelics and left-aligns upstream, so each record
    carries a single ALT and the comparison is unambiguous:
      - SNV: equal length 1
      - MNV: equal length > 1
      - INS: ALT longer than REF
      - DEL: REF longer than ALT
    """
    len_ref, len_alt = len(ref), len(alt)
    if len_ref == len_alt:
        return "SNV" if len_ref == 1 else "MNV"
    return "INS" if len_alt > len_ref else "DEL"


def _scalar_int(value):
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        return None
    return int(value)


def _ad_pair(ad):
    if ad is None:
        return None, None
    if not isinstance(ad, (tuple, list)):
        return None, None
    ref = int(ad[0]) if len(ad) > 0 and ad[0] is not None else None
    alt = int(ad[1]) if len(ad) > 1 and ad[1] is not None else None
    return ref, alt


def _scalar_float(value):
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        if not value or value[0] is None:
            return None
        return float(value[0])
    return float(value)


def main():
    vcf_path = sys.argv[1]
    sample_name = sys.argv[2]
    out_path = sys.argv[3]

    rows = []
    vcf = pysam.VariantFile(vcf_path)
    callers = list(vcf.header.samples)

    for rec in vcf:
        # BCFTOOLS_NORM split multi-allelics upstream, so ALT is single.
        alt = rec.alts[0] if rec.alts else None
        if alt is None:
            continue

        for caller in callers:
            sample = rec.samples[caller]
            dp = _scalar_int(sample.get("DP"))
            if dp is None or dp == 0:
                continue
            ad_ref, ad_alt = _ad_pair(sample.get("AD"))
            af = _scalar_float(sample.get("VAF"))
            caller_qual = _scalar_float(sample.get("VARQUAL"))

            rows.append(
                {
                    "sample_name": sample_name,
                    "caller": caller,
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": alt,
                    "type": _variant_type(rec.ref, alt),
                    "qual": caller_qual,
                    "dp": dp,
                    "ad_ref": ad_ref,
                    "ad_alt": ad_alt,
                    "af": af,
                    "caller_specific": "{}",
                }
            )

    df = pl.DataFrame(rows, schema=SCHEMA)
    df.write_parquet(out_path)


if __name__ == "__main__":
    main()
