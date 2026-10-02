#!/usr/bin/env python3
"""Migrate the plain amplicon FASTA to a self-describing enriched FASTA.

Reads the panel design CSV and the plain alignment FASTA.
Writes an enriched FASTA whose genomic metadata lives in the record description.
That is after the FASTA ID; everything after the first space on the `>` line.
The record ID will be left unchanged. This ID will be join key for example at the @SQ
contig in the BAM, the `feature` key in the downstream analysis, for example count_variant.py
as well as what input_validation.nf validates as `[A-Z_0-9]+`).

Aligners truncate the contig name at the first whitespace, so the description will be
invisible to alignment. One verification that can be is to check that the enriched FASTA
must yield identical @SQ / identical alignments.

Header format of the output FASTA is shown in the example below:

    >RS10750803 legacy_name=RS10750803_1_63 chrom=chr11 g0=67734055 amplicon_start=67734156 amplicon_end=67734512 snp_pos=67734444

The contig ID will be shortend to extract the first token after split by "_"
for example `RS10750803`. The legacy `_1_<n>` index form is dropped from the ID
and preserved in the description as `legacy_name=`. This is intentional design
choice to make the ID more traceible and being able to backtrack.

The short ID still satisfies input_validation's `[A-Z_0-9]+`, and `contig_to_rsid` still resolves it.

Stored fields, everything else derived (lossless, no stored redundancy):

- legacy_name=    the original `RS..._1_<n>` contig name (traceability only)
- chrom=          real chromosome
- g0=             1-based genomic coord aligned to amplicon-local 0-based pos 0
                  (== TSV contig_start == Amplicon_Start - LEFT_PAD == offset + 1).
                  The single lift anchor:  gpos = g0 + pos_local_0based  (NO +1).
- amplicon_start= targeted insert start (Amplicon_Start)
- amplicon_end=   targeted insert end   (Amplicon_End)
- snp_pos=        designed SNP genomic position (on-target classification); OPTIONAL,
                  omitted from the header when the design CSV has no snps_pos for the amplicon

Note for the developers: `fasta_len` is deliberately NOT stored; it is len(seq)
at read time and would drift from the sequence. `rsid`, `offset`, `contig_end`,
`Amplicon_Length` are all derivable and likewise not stored.

This generator intentionally uses only the stdlib and NOT biopython:
FASTA parsing here is trivial, it keeps the off-pipeline tool dependency-free.
This is because of Wave container in the pipeline doesn't have biopython installed.
So this will parse `>` lines directly. The locked thing is the header *format*,
not the parser library.

Output is written single-line-per-sequence to match the existing FASTA layout, so the
only diff vs the current file is the header lines (makes the alignment-identity argument
trivially easy).

Usage:
    make_self_describing_fasta.py --csv CSV --fasta PLAIN_FASTA --out ENRICHED_FASTA
                                  [--left-pad 101]
"""

import argparse
import csv
import logging
import re
import sys

logger = logging.getLogger("make_self_describing_fasta")

# input_validation.nf validates the @SQ contig name against this pattern.
# so we need to have it here as well.
ID_RE = re.compile(r"^[A-Z_0-9]+$")


def read_amplicons(csv_path):
    """rsid -> dict(chrom,start,end,length,snps)."""
    amp = {}
    with open(csv_path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            snp_raw = (r.get("snps_pos") or "").strip()
            amp[r["Name"].lower()] = dict(
                chrom=r["Chrom"],
                start=int(r["Amplicon_Start"]),
                end=int(r["Amplicon_End"]),
                length=int(r["Amplicon_Length"]),
                snps=int(snp_raw) if snp_raw else None,
            )
    return amp


def contig_to_rsid(contig):
    """`RS1604513_1_49` -> `rs1604513` (token before the first underscore)."""
    return contig.split("_")[0].lower()


def read_fasta_records(path):
    """Yield (contig_id, sequence) in file order."""
    contig, chunks = None, []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if contig is not None:
                    yield contig, "".join(chunks)
                contig = line[1:].split()[0]  # ID = token up to first whitespace
                chunks = []
            elif line:
                chunks.append(line.strip())
    if contig is not None:
        yield contig, "".join(chunks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="panel design CSV")
    ap.add_argument("--fasta", required=True, help="plain alignment FASTA (bare IDs)")
    ap.add_argument("--out", required=True, help="enriched FASTA to write")
    ap.add_argument("--left-pad", type=int, default=101)
    a = ap.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    amp = read_amplicons(a.csv)

    n = 0
    bad = []
    # store new short ID to legacy contig mapping that produced it
    seen_ids = {}
    with open(a.out, "w") as out:
        # this is not efficient obviously, but works for a small FASTA file
        # for amplicons and limited number of records in the reference it should
        # do just fine!
        for contig, seq in read_fasta_records(a.fasta):
            rs = contig_to_rsid(contig)
            if rs not in amp:
                raise KeyError(f"contig {contig} has no CSV row {rs}")
            row = amp[rs]
            start, end, length, snp = (
                row["start"],
                row["end"],
                row["length"],
                row["snps"],
            )
            g0 = start - a.left_pad

            # split happens here for the id
            # RS10750803_1_63 -> RS10750803
            new_id = contig.split("_")[0]

            # bunch of conditions for:
            # - defenses against validations
            # - id collisions
            # - amplicon length mismatches
            # - snp pos being outside the ROI
            errs = []
            if not ID_RE.match(new_id):
                errs.append(f"short ID {new_id!r} fails input_validation [A-Z_0-9]+")
            if new_id in seen_ids:
                errs.append(f"short ID {new_id} collides with {seen_ids[new_id]}")
            if end - start + 1 != length:
                errs.append(
                    f"CSV Amplicon_Length {length} != end-start+1 {end - start + 1}"
                )
            if len(seq) != length + 201:
                errs.append(
                    f"len(seq) {len(seq)} != Amplicon_Length+201 {length + 201}"
                )
            if snp is not None and not (start <= snp <= end):
                errs.append(f"snp_pos {snp} outside insert [{start},{end}]")
            if errs:
                bad.append((contig, errs))
            seen_ids[new_id] = contig

            desc = (
                f"legacy_name={contig} chrom={row['chrom']} g0={g0} "
                f"amplicon_start={start} amplicon_end={end}"
            )
            if snp is not None:
                desc += f" snp_pos={snp}"
            out.write(f">{new_id} {desc}\n{seq}\n")
            n += 1

    if bad:
        logger.error("FAILED: %d contig(s) violated a generation invariant:", len(bad))
        for contig, errs in bad:
            for e in errs:
                logger.error("    %s: %s", contig, e)
        logger.error("enriched FASTA was written but MUST NOT be used; fix inputs.")
        return 1

    logger.info("wrote %d enriched records -> %s", n, a.out)
    logger.info("all generation invariants passed (pad, length, snp-in-insert)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
