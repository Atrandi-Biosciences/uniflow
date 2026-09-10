# Variant caller merge (COMB-439)

## Goal

Combine per-caller VCFs (vardict, freebayes, bcftools, gatkhc) for a sample
into one consensus VCF with cross-caller INFO summaries and per-caller
FORMAT columns. SNV-only. No quality filtering in this PR — that lands as
a separate upstream module later.

Lofreq is wired through but excluded from the default `variant_callers`:
its site-only output has no `FORMAT/AD,DP`, so it cannot contribute to the
consensus FORMAT-derived summaries without a custom INFO→FORMAT lift. The
template + merge plumbing remains in the codebase for future revival; the
lift itself is parked because no pure-bcftools path exists (`+fill-tags`
cannot construct multi-value FORMAT fields from INFO) and the available
alternatives (pysam script, awk text rewrite) carry more custom code than
the immediate value justifies.

## Strategy vs sarek 3.8.1

Sarek keeps per-caller VCFs separate, uses `INFO/SOURCE=<filename>` for
provenance, and does not harmonize VAF. We diverge: we actually merge into
one VCF per sample with caller-as-pseudo-sample columns, and we expose the
across-caller VAF spread as `INFO/VAF_{min,max}`. We adopt sarek's
`conf/modules.config` + `ext.args` pattern and its
`bcftools norm --multiallelics -both --rm-dup all` flags.

## Decisions

- **Merge policy:** union. Caller identity lives in the per-sample columns
  (sample name = caller string), so no `INFO/CALLERS` string list is needed —
  downstream consumers derive caller membership from `FORMAT/DP>0` per
  sample. `INFO/NCALLERS` is kept as a precomputed cross-caller count.
- **Consensus INFO surface:** `INFO/VAF_{min,max}`, `INFO/NCALLERS`,
  `INFO/DP_{min,median,max}`, `INFO/AD_REF_{min,median,max}`,
  `INFO/AD_ALT_{min,median,max}`. All derived post-merge from the per-caller
  FORMAT columns via a bcftools-native `+fill-tags` pipe — no awk, no custom
  scripting. `FORMAT` is restricted to `VAF / AD / DP / GT` so the consensus
  VCF lands on a fixed schema for tabular conversion (CSV/parquet/H5)
  downstream. Counts (DP, AD) report min/median/max; fractional VAF reports
  only min/max (bcftools median over an even count returns a float that
  doesn't round-trip cleanly into a Number=1,Type=Integer field; for VAF the
  spread is the load-bearing signal anyway).
- **VAF exposure:** `FORMAT/VAF` is computed post-merge from `FORMAT/AD`
  (REF, ALT) by `+fill-tags`'s built-in VAF rule, so each per-caller sample
  column carries its own VAF. Across-caller spread is exposed as
  `INFO/VAF_{min,max}`. No per-caller VAF dispatch before merge.
- **Sample column naming:** caller name only (`lofreq`, `vardict`, ...).
  Sample identity lives in filename + `meta.sample_name`.
- **Step order:** norm → fill-VARQUAL → (lofreq only: attach sample column)
  → reheader sample → merge → annotate consensus (FORMAT/VAF and all INFO
  summaries computed here).
- **Config pattern:** introduce `conf/modules.config` with `withName` +
  `ext.args`. Only new/touched modules adopt it; existing modules stay as-is.
- **QC step:** ship in this PR. Nextflow port of
  `docs/adr/lofreq/compare_caller_runs.sh`, run per sample.

## VARQUAL: preserving per-caller site quality through merge

`bcftools merge` collapses each caller's QUAL into a single merged QUAL,
losing the per-caller value. Two changes together preserve per-caller QUAL
through merge:

1. `BCFTOOLS_FILL_TAGS` fills `INFO/VARQUAL` with the caller's QUAL but
   declares the field with Number = the count of active callers
   (`-t 'INFO/VARQUAL:${n_callers}=QUAL'`, where `n_callers` is derived at
   runtime from `params.variant_callers.tokenize(',').size()`). The
   per-caller VCF carries a single value under a length-`n_callers` header
   declaration. Subsetting / adding callers updates the declaration
   automatically.
2. `BCFTOOLS_MERGE` is invoked with `--info-rules VARQUAL:join` (set in
   `conf/modules.config`). The `join` rule concatenates each input file's
   `INFO/VARQUAL` into a per-record list in input-file order — for the
   per-sample merge that yields a length-`n_callers` list aligned with the
   sample columns, with `.` in the slot of any caller missing the site.
   Without `--info-rules` the default behaviour is first-non-missing
   (lossy).

`BCFTOOLS_ANNOTATE_CONSENSUS` keeps `INFO/VARQUAL` through the consensus
build via `bcftools annotate --remove '^INFO/VAF,^INFO/VARQUAL'` on the
leading strip.

VAF is *not* filled per-caller — `BCFTOOLS_ANNOTATE_CONSENSUS` computes
`FORMAT/VAF` from `FORMAT/AD` post-merge using `bcftools +fill-tags`'s
built-in rule, avoiding the need for per-caller VAF dispatch.

## Lofreq sample-column attachment (dormant)

Lofreq is dormant in the default config but its plumbing remains wired so
the path can be reactivated with one config flip.

Lofreq's output is site-only (`#CHROM POS ID REF ALT QUAL FILTER INFO` — no
`FORMAT` block, no sample column). `INFO/VARQUAL` can be filled directly
because the expression is QUAL→INFO, but the downstream
`BCFTOOLS_REHEADER_SAMPLE` refuses to rename a non-existent sample column
("missing FORMAT fields, cowardly refusing to add samples"), so we attach
one before reheader.

The mechanism: `TEMPLATE_VCF` emits a header-only VCF declaring `FORMAT/GT`
and a placeholder sample column named after `meta.sample_name`. A second
`bcftools merge` invocation (`BCFTOOLS_MERGE_LOFREQ`, aliased from
`BCFTOOLS_MERGE`) merges the lofreq VCF with the template, producing
records that inherit lofreq's `INFO` (including `INFO/VARQUAL`) plus a
`FORMAT/GT=./.` sample column. The placeholder sample name is then
overwritten by `BCFTOOLS_REHEADER_SAMPLE` to `lofreq` like every other
caller.

What's missing for full participation: `FORMAT/AD` and `FORMAT/DP` are not
populated, so even when lofreq is re-enabled it contributes only `GT=./.`
and its `INFO/*` to the merge — not the FORMAT-derived consensus summaries
(`INFO/DP_*`, `INFO/AD_*`, `FORMAT/VAF`). Closing this gap requires
deriving `FORMAT/DP` from `INFO/DP` and `FORMAT/AD` from `INFO/DP4` per
record. `bcftools +fill-tags` cannot construct multi-value FORMAT fields
from same-record INFO (verified empirically); `vembrane` is for filter /
external annotation / table output, not field derivation. A targeted
pysam (or awk) script would do it cleanly but is deferred until lofreq's
value justifies the custom code.

## Modules

New under `modules/`:

- `bcftools/fill_tags/` — `INFO/VARQUAL:1=QUAL` (uniform across callers);
  preserves each caller's site QUAL through the merge
- `bcftools/reheader_sample/` — rename single sample to caller string
- `bcftools/merge/` — generalised `bcftools merge --force-samples -m none`
  with `ext.prefix` override; used both for per-caller merge in
  `variant_merging.nf` and for lofreq+template merge in `variant_harmonize.nf`
  (via the `BCFTOOLS_MERGE_LOFREQ` alias)
- `bcftools/annotate_consensus/` — strips all per-caller INFO from the
  merge and rebuilds the consensus surface via chained `+fill-tags`:
  `INFO/VAF_{min,max}`, `INFO/NCALLERS`, `INFO/DP_{min,median,max}`,
  `INFO/AD_REF_{min,median,max}`, `INFO/AD_ALT_{min,median,max}`. Closes
  by restricting FORMAT to `VAF / AD / DP / GT`. Pinned to bcftools 1.23.1
  for the `FORMAT/AD[:0,]` / `FORMAT/AD[:1,]` subscript syntax (same
  constraint as `bcftools/fill_tags/`).
- `template_vcf/` — header-only VCF generator (sample column with
  `FORMAT/GT`), used by the lofreq merge path
- `qc/caller_compare/` — Nextflow wrap of compare_caller_runs.sh

Edited:

- `bcftools/norm/main.nf` — switch hardcoded flags to `${task.ext.args}`;
  conf supplies `--multiallelics -both --rm-dup all`. SNV-only `bcftools view`
  pipe stays (pipeline policy, not configurable).

Moved/renamed:

- `modules/lofreq/main.nf` → `modules/lofreq/call/main.nf`; process
  `LOFREQ` → `LOFREQ_CALL`. Aligns lofreq's layout with the per-subcommand
  convention used by other multi-step tool families.

## Pipeline

For each per-caller VCF: `BCFTOOLS_NORM` splits multiallelics and dedups,
then `BCFTOOLS_FILL_TAGS` copies `QUAL` into `INFO/VARQUAL`. Lofreq's
site-only VCF then passes through `TEMPLATE_VCF` + `BCFTOOLS_MERGE_LOFREQ`
to gain a sample column; the other four callers already have one.

The harmonised VCFs feed two consumers, grouped by `meta.sample_name`:
`CALLER_COMPARE_QC` produces the per-sample comparison report, and the merge
path continues through `BCFTOOLS_REHEADER_SAMPLE` (rename the single sample
column to the caller string) into `BCFTOOLS_MERGE`. The merged VCF then
passes through `BCFTOOLS_ANNOTATE_CONSENSUS`, which strips the per-caller
INFO carried over by merge and synthesises the consensus surface
(`INFO/VAF_{min,max}`, `INFO/NCALLERS`, `INFO/DP_{min,median,max}`,
`INFO/AD_REF_{min,median,max}`, `INFO/AD_ALT_{min,median,max}`) from the
per-caller `FORMAT` columns via chained `bcftools +fill-tags` calls.
`FORMAT` is restricted to `VAF / AD / DP / GT` on the consensus output.

## Out of scope

- Quality filtering (next module; will sit upstream of merge).
- Indel handling — pipeline is SNV-only.
- Cross-sample merge (single sample per output VCF).
- Lofreq INFO→FORMAT lift (deferred — lofreq dormant in the default
  `variant_callers`; modules retained for future revival).
