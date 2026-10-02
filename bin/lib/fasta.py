"""Shared FASTA parser and utilities."""

import pysam


def parse_fasta(path):
    """Yield (record_id, fields, bad_tokens, sequence) for each record.

    - The description is everything after the first whitespace, split into key=value pairs.
    Keys are lowercased so AS=/as= and M5=/m5= are the same field. Tokens without an
    = are returned separately rather than dropped, so callers can reject a malformed header
    instead of silently ignoring it.

    - Lowercased fields are considered valid so AS=/as= and M5=/m5=

    record_id is the first whitespace-delimited token after > (may be "" for a bare
    > line), fields maps lowercased key to value, bad_tokens lists description
    tokens that carried no =, and sequence is the concatenated sequence lines ("" for
    a record that has none).
    """
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            fields, bad_tokens = {}, []
            for token in (entry.comment or "").split():
                if "=" in token:
                    key, value = token.split("=", 1)
                    fields[key.lower()] = value
                else:
                    bad_tokens.append(token)
            yield entry.name, fields, bad_tokens, entry.sequence
