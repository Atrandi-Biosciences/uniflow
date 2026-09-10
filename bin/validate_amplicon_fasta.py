#!/usr/bin/env python3
"""Validate an amplicon FASTA against the amplicon FASTA specification.

Panel-agnostic runtime check. Every description field is optional, so a plain FASTA with no
`key=value` description will validate. Whatever is declared is checked for being well formed and
self-consistent according to the RFC.

Errors exit nonzero and fail the pipeline run, whereas warnings are reported and the run
continues:

  errors    a description token that is not `key=value`; an ID that is missing or falls
            outside [A-Za-z0-9_]; a duplicate ID; a record with no sequence; a base outside
            ACGTN; a g0/amplicon_start/amplicon_end/snp_pos that is not a non-negative
            integer; amplicon_end below amplicon_start; snp_pos outside the amplicon; an M5
            that does not match the MD5 of the sequence

  warnings  g0 past amplicon_start (a padded contig cannot contain its own insert);
            lowercase (soft-masked) bases; an unrecognised description key; two records
            sharing at least SIMILARITY_WARN of their k-mers, an ambiguity the panel
            cannot resolve downstream and which splits their per-cell counts

Cross-field checks run only when both fields are present. The whole file is checked before
exiting, so one run will report every problem rather than only the first.

Besides the report on stderr, every amplicon pair is streamed to
`amplicon_fasta_similarity.csv` (ref_1, ref_2, similarity) in the working directory, warning or
not, so the pairs below the warning threshold stay visible instead of being discarded.

Usage:
    validate_amplicon_fasta.py FASTA
"""

import argparse
import hashlib
import logging
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations

from lib.fasta import parse_fasta

logger = logging.getLogger("validate_amplicon_fasta")

ERROR = "error"
WARNING = "warning"

VALID_ID = re.compile(r"^[A-Za-z0-9_]+$")
KNOWN_KEYS = frozenset(
    {
        "as",
        "chrom",
        "g0",
        "amplicon_start",
        "amplicon_end",
        "snp_pos",
        "legacy_name",
        "m5",
    }
)
VALID_BASES = frozenset("ACGTN")
COORD_KEYS = ("g0", "amplicon_start", "amplicon_end")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# A 31-mer is long enough that a shared hit means real sequence identity rather than
# chance. Measured across 208 amplicons in our own panels, the worst unrelated pair
# shares 1.9% of its 31-mers, and three of the four panels share none at all.
#
# The threshold is a property of the panel.
# What makes a pair of amplicons references dangerous is a shared stretch long enough
# to swallow a whole insert, since only then is there nothing unique left to place
# that insert by. If we assume that the containment for a shared block is roughly
# block_length / amplicon_length, so a 150bp stretch of a 400bp amplicon reference
# is ~0.34.
#
# So a similarity warning of 0.3 is just under that also catches scattered divergences!
# The 150bp is an assumed property of the assay, not an input of course.
# If inserts become shorter the threshold has to be adjusted with them
KMER_SIZE = 31
SIMILARITY_WARN = 0.3

# Written beside the log, in the process working directory
SIMILARITY_CSV = "amplicon_fasta_similarity.csv"


def canonical_kmers(seq):
    """Generate strand-agnostic k-mer set

    each k-mer stored as min(itself, its reverse complement)."""
    kmers = set()
    for i in range(len(seq) - KMER_SIZE + 1):
        kmer = seq[i : i + KMER_SIZE]
        kmers.add(min(kmer, kmer.translate(COMPLEMENT)[::-1]))
    return kmers


def check_record(record_id, fields, bad_tokens, seq):
    """Yield (level, record_id, message) for every single-record rule that fires."""
    for token in bad_tokens:
        yield ERROR, record_id, f"description token {token!r} is not key=value"
    if not VALID_ID.match(record_id):
        yield ERROR, record_id, f"ID {record_id!r} is not alphanumeric or underscore"
    for key in sorted(set(fields) - KNOWN_KEYS):
        yield WARNING, record_id, f"unrecognised description key {key!r}"

    if not seq:
        yield ERROR, record_id, "record has no sequence"

    upper = seq.upper()
    for base in sorted(set(upper) - VALID_BASES):
        yield ERROR, record_id, f"invalid base {base!r} in sequence"
    if seq != upper:
        yield WARNING, record_id, "sequence contains lowercase (soft-masked) bases"

    coords = {}
    for key in COORD_KEYS:
        if key in fields:
            value = fields[key]
            if value.isdecimal():
                coords[key] = int(value)
            else:
                yield ERROR, record_id, f"{key}={value!r} is not a non-negative integer"
    snps = []
    for value in fields.get("snp_pos", "").split("|"):
        if not value:
            continue
        if value.isdecimal():
            snps.append(int(value))
        else:
            yield ERROR, record_id, (
                f"snp_pos element {value!r} is not a non-negative integer"
            )

    start, end = coords.get("amplicon_start"), coords.get("amplicon_end")
    if start is not None and end is not None and end < start:
        yield ERROR, record_id, f"amplicon_end {end} < amplicon_start {start}"
    elif start is not None and end is not None:
        for snp in snps:
            if not start <= snp <= end:
                yield ERROR, record_id, f"snp_pos {snp} outside amplicon [{start},{end}]"
    if "g0" in coords and start is not None and coords["g0"] > start:
        yield WARNING, record_id, f"g0 {coords['g0']} > amplicon_start {start}"

    if "m5" in fields:
        digest = hashlib.md5(upper.encode()).hexdigest()
        if digest != fields["m5"].lower():
            yield ERROR, record_id, f"M5 {fields['m5']} != sequence MD5 {digest}"


def check_similarity(records, out):
    """Check similarity of pairs of sequences in the input FASTA

    Note: This will write every pair's containment to `out` and yield a warning for those
    at the threshold mentioned above.

    A row is written as each pair is scored, so the table records what the warning will
    hide: the pairs below SIMILARITY_WARN, which are the ones that say whether the
    threshold sits above the panel's lower noise level or in the middle of it. Pairs sharing
    no k-mer at all are not scored and have no row; their similarity is 0 by construction.

    Note 2: Not writing all the FASTA similarity is just to make this code simpler
    """
    out.write("ref_1,ref_2,similarity\n")
    kmer_sets = [canonical_kmers(seq.upper()) for _, seq in records]
    index = defaultdict(list)
    for i, kmers in enumerate(kmer_sets):
        for kmer in kmers:
            index[kmer].append(i)

    shared = Counter()
    for owners in index.values():
        if len(owners) > 1:
            shared.update(combinations(owners, 2))

    for (i, j), count in sorted(shared.items()):
        # Containment is over the smaller set. Swap i,j
        if len(kmer_sets[j]) < len(kmer_sets[i]):
            i, j = j, i
        containment = count / len(kmer_sets[i])
        out.write(f"{records[i][0]},{records[j][0]},{containment:.6f}\n")
        if containment >= SIMILARITY_WARN:
            yield WARNING, records[i][0], (
                f"shares {containment:.1%} of its {KMER_SIZE}-mers with "
                f"{records[j][0]}; reads _might_ align ambiguously between them"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fasta", help="amplicon FASTA to validate")
    a = ap.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    # Collect everything before reporting so one run generates every problem in
    # the input FASTA
    problems, records = [], []
    for record_id, fields, bad_tokens, seq in parse_fasta(a.fasta):
        problems.extend(check_record(record_id, fields, bad_tokens, seq))
        records.append((record_id, seq))

    if not records:
        problems.append((ERROR, "-", "file contains no FASTA records"))

    for record_id, count in sorted(Counter(r[0] for r in records).items()):
        if count > 1:
            problems.append((ERROR, record_id, f"duplicate ID, appears {count} times"))

    with open(SIMILARITY_CSV, "w") as out:
        problems.extend(check_similarity([r for r in records if r[1]], out))

    errors = [(rid, message) for level, rid, message in problems if level == ERROR]
    warnings = [(rid, message) for level, rid, message in problems if level == WARNING]

    logger.info("records: %d", len(records))
    for record_id, message in warnings:
        logger.warning("    %s: %s", record_id, message)
    if errors:
        logger.error("FAILED: %d problem(s):", len(errors))
        for record_id, message in errors:
            logger.error("    %s: %s", record_id, message)
        return 1

    logger.info("PASS: %d records, %d warning(s)", len(records), len(warnings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
