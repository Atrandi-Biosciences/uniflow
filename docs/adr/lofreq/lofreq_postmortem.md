# LoFreq Postmortem

| Field | Value |
|---|---|
| Status | Paused. LoFreq is excluded from the default caller set; the implemented module remains in-tree at [modules/lofreq/main.nf](../../modules/lofreq/main.nf) as the entry point if revived |
| Author | Hassan Foroughi Asl |
| Last updated | 2026-05-15 |
| Related code | [modules/lofreq/main.nf](../../modules/lofreq/main.nf), merged via `feature/COMB-416-implement-call_lofreq-module` (commit a0301b7, PR #71) |
| Related docs | [docs/variant_calling_specs.md](../variant_calling_specs.md) ;  overall variant-calling design |

> **Path note (2026-05-19):** module paths in this document refer to the
> codebase as of the postmortem date. The LoFreq module has since been moved
> from `modules/lofreq/main.nf` (`LOFREQ`) to `modules/lofreq/call/main.nf`
> (`LOFREQ_CALL`) to match the per-subcommand layout used by other
> multi-step tool families. Links below may resolve to the old path.

## TL;DR

LoFreq's exact Poisson-binomial DP is `O(N × K)` per locus, where `N` is per-locus coverage and `K` is `max(noncons_counts)`. noncons_counts refers to the bases that are non-reference. At our bulk-aggregated per-amplicon depth (≈ 650K read depth), loci with a real high-VAF variant hit ~10¹¹ tests in a for loop each and this takes hours. There is no algorithmic pruning as far as I could tell or it is probably asymmetric; meaning it will only work well on non-variant locus.

LoFreq has `-t/--approx-threshold`, which switches to a Poisson approximation above a depth threshold and collapses the cost to `O(N)` (according to their paper and source code). It is off by default. With `-t 10000` plus a selection of option set, LoFreq is usable on this data. I suggest nevertheless excluding LoFreq or at least having it as an optional variant caller, because bcftools / freebayes / vardict / gatk-hc already cover germline SNV calling without this problem; the option set below is the entry point if we change this decision in the future.

I ran a separate diagnostic to check the BAM file produced by Minimap2, documented below; but this doesn't eliviate the issue with LoFreq's performance.

## 1. Context

### 1.1 Pipeline

The variant-calling subworkflow (see [docs/variant_calling_specs.md](../variant_calling_specs.md)) runs multiple germline SNV callers in parallel on amplicon BAMs and produces per-caller normalised VCFs. The candidate caller set under evaluation is: bcftools, freebayes, vardict, gatk-hc, **lofreq**. This document covers only the LoFreq.

### 1.2 Sample and chemistry

- **Test sample:** RAR143 S11, downloaded from Latch via `latch cp`.
- **Chemistry:** single-cell amplicon. R1 carries the genomic insert; **R2 carries only the cell barcode** and is not aligned to the reference. The aligned BAM is therefore single-end-shaped (`RNEXT=*`, `TLEN=0` everywhere) ;  this is by design. Per-amplicon depths in this document are the bulk aggregate across cells; per-cell depth at any position is not evaluated.
- **Aligner:** minimap2, `-ax sr` preset, amplicon-as-contig reference (each amplicon is its own short contig, e.g. `RS5747378_1_57`).

### 1.3 Software versions

- LoFreq: 2.1.5 ;  `quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4` (pinned in [modules/lofreq/main.nf](../../modules/lofreq/main.nf)).
- Reference upstream: [github.com/CSB5/lofreq](https://github.com/CSB5/lofreq) ;  `src/lofreq/snpcaller.c`, `src/lofreq/lofreq_call.c`, `src/lofreq/defaults.h`. Where the Poisson approximation and the methodology are described.

### 1.4 Command that exposed the issue

```
lofreq call-parallel --no-ext-baq --force-overwrite --debug --pp-threads 6 \
  --bed intergenic_amplicon_50.bed --verbose \
  -f intergenic_amplicon_50.fasta \
  -o amplicon-panel__protK_yes_1_min_ext_rep1.lofreq.vcf.gz \
  aligned.sorted.bam
```

Run time on the full data was unusable. It stalled for 8 hours on Latchbio, and it had to be manually aborted. I investigated individual options and a way to fine tune and speed it up. See § 5 below.

## 2. Symptom

Three consecutive positions on amplicon `RS5747378_1_57` from a `lofreq call --debug` run:

```
DEBUG(plp.c|compile_plp_col): Processing RS5747378_1_57:347
RS5747378_1_57	347	T	T	counts:rv/fw A:0/4 C:0/168 G:0/1364 T:0/602116 N:0/0 heads:0 tails:1 ins:0 del:6 hrun=1
DEBUG(lofreq_call.c|call_snvs): RS5747378_1_57 347: passing down 603652 quals with noncons_counts (4, 168, 1364) to snpcaller(num_snv_tests=327 conf->bonf=327, conf->sig=0.010000)

DEBUG(plp.c|compile_plp_col): Processing RS5747378_1_57:348
RS5747378_1_57	348	C	G	counts:rv/fw A:0/187 C:0/170711 G:0/432487 T:0/268 N:0/0 heads:1 tails:0 ins:1 del:5 hrun=5
DEBUG(lofreq_call.c|call_snvs): RS5747378_1_57 348: passing down 603653 quals with noncons_counts (187, 432487, 268) to snpcaller(num_snv_tests=330 conf->bonf=330, conf->sig=0.010000)
DEBUG(lofreq_call.c|call_snvs): low freq snp: RS5747378_1_57 348 G>G pv-prob:3.3621e-4932;pv-qual:49314 counts-raw:432487/603663=0.716438 counts-filt:432487/603653=0.716450

DEBUG(plp.c|compile_plp_col): Processing RS5747378_1_57:349
RS5747378_1_57	349	A	G	counts:rv/fw A:0/170965 C:0/10 G:0/432660 T:0/18 N:0/0 heads:5 tails:0 ins:1 del:0 hrun=5
DEBUG(lofreq_call.c|call_snvs): RS5747378_1_57 349: passing down 603653 quals with noncons_counts (10, 432660, 18) to snpcaller(num_snv_tests=333 conf->bonf=333, conf->sig=0.010000)
```

Reading the fields:

- `counts:rv/fw A:.. C:.. G:.. T:..` ;  strand-split per-base read counts.
- `noncons_counts (x, y, z)` ;  the three non-reference base counts (reference base dropped from A/C/G/T, slot order shifts with the reference):
  - Pos 347 (ref `T`): `(A, C, G) = (4, 168, 1364)` ;  non-variant locus, noise only.
  - Pos 348 (ref `C`): `(A, G, T) = (187, 432487, 268)` ;  apparent G variant at 71.6% VAF.
  - Pos 349 (ref `A`): `(C, G, T) = (10, 432660, 18)` ;  apparent G variant at 71.7% VAF.
- `num_snv_tests` / `conf->bonf` ;  dynamic Bonferroni factor, +3 per informative locus (327, 330, and 333 respectively).
- `pv-prob:3.3621e-4932` ;  final p-value at pos 348. After multiplying by `bonf=330` it is still ~10⁻⁴⁹²⁹, ~4,927 orders of magnitude below `sig=0.01`. The pruning check `pvalue * bonf > sig` is never satisfied, thus there will be no early exit, and it will cause a the full `O(N × K)` DP runs.

## 3. Root cause

LoFreq computes a per-locus p-value under the null "all `K = max_noncons_count` non-reference bases are sequencing errors" using an exact Poisson-binomial dynamic program. The DP lives in `pruned_calc_prob_dist` in [src/lofreq/snpcaller.c](https://github.com/CSB5/lofreq/blob/master/src/lofreq/snpcaller.c) and is wrapped by `poissbin`, called from `snpcaller`.

```c
for (n=1; n<=N; n++) { // N = num_err_probs = locus coverage
    for (k=MIN(n, K-1); k>=1; k--) { // K = max_noncons_count
        probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                             probvec_prev[k-1] + log_pn);
    }
}
```

According to the comments in the code and the publication:
> the total inner-loop count is `Σ min(n, K-1) ≈ K × N` for `K ≪ N`, with a triangular correction when `K ~ N`
thus the the early-exit will be

```c
if (pvalue * (double)bonf_factor > sig_level) return probvec;
```

The running p-value (probability of observing ≥ K errors or mutations in the first `n` reads) starts tiny and grows toward the final value. The exit only happens when the running p-value *grows past* `sig / bonf`. For non-variant loci (small `K`, p-value climbs past threshold quickly) this works and is fast. For real high-VAF variants the p-value stays astronomically small as seen above(`3.36e-4932`), the threshold is never crossed, and the loop runs to completion; thus it takes forever.

Observed on `RS5747378_1_57`:

| Pos | N       | K       | Inner-loop ops (`Σ min(n, K-1)`) | Notes              |
|-----|---------|---------|----------------------------------|--------------------|
| 347 | 603,652 | 1,364   | ~8.2 × 10⁸                       | noise              |
| 348 | 603,653 | 432,487 | ~1.7 × 10¹¹                      | real G variant     |
| 349 | 603,653 | 432,660 | ~1.7 × 10¹¹                      | real G variant     |

Two loci alone account for ~3.4 × 10¹¹ scalar `log_sum` ops.

Hardware does not help meaningfully:

- `lofreq call-parallel` splits by region. A single high VAF locus still runs single-threaded.
- The inner loop is non-vectorisable: `log_sum(a, b) = max(a, b) + log1p(exp(-|a-b|))` is branchy, and each `k` depends on `k-1` of the previous `n`.
- Best realistic CPU/SIMD speedup is ~2×
- The `-t` algorithmic fix is ~10⁵× (see below § 5)

## 4. MAPQ side quest

I took a detour to investigate whether the slow amplicon was due to the mapping issues, in case we could use `-s/--src-qual` option in `lofreq` . It turned out not to be.

### 4.1 MAPQ distribution

```bash
samtools view aligned.sorted.bam | awk '{print $5}' | sort -n | uniq -c
```

Tail (11,111,851 reads total):

```
   2902 49
   3611 50
   9257 51
 225881 52
  88353 53
   1944 54
   1854 55
   2031 56
   2234 57
   2216 58
   7940 59
10216037 60
```

92% of reads at MAPQ 60. Distinct spike at MAPQ 52. This is what justifies `--min-mq 50` (keeps the 52/53 cluster; drops only the ambiguous tail).

### 4.2 MAPQ-52 spike

MAPQ-52 concentrates on one amplicon, which is coincidentally the region that LoFreq stalls for hours:

``` bash
samtools view -q 52 aligned.sorted.bam | awk '$5==52 {print $3}' | sort | uniq -c | sort -rn | head

 223348 RS5747378_1_57
    506 RS12147394_1_26
    394 RS9786182_1_46
    ...
```

99% of the MAPQ-52 reads are on the same amplicon that hangs LoFreq. Initial hypothesis: paralog-driven competing alignment? Difficult to map region?

### 4.3 Tag inspection rules out paralog

```bash
samtools view -q 52 aligned.sorted.bam RS5747378_1_57 | head -1 | tr '\t' '\n' | tail -n +6
```

Representative read:

```bash
148M ... NM:i:1  ms:i:286  AS:i:286  nn:i:0  tp:A:P  cm:i:20
        s1:i:127  s2:i:0  de:f:0.0068  rl:i:0
```

Confirmed across all reads on this amplicon:

```bash
samtools view -q 52 aligned.sorted.bam RS5747378_1_57 | tr '\t' '\n' | grep 's2:i:' | grep -v "s2:i:0" | head
# (returns nothing)
```

`s2:i:0` everywhere ; minimap2 found **no competing chain**. Clean `148M` alignment, single mismatch, no Ns in reference, no soft-clipping. The MAPQ cap at 52 comes from minimap2's `fac` factor in its MAPQ formula (chain confidence relative to alignment length on a short amplicon-as-contig). It is probably artifact of the panel design, not a paralog signal.

Consequences:

- `-s/--src-qual` is not useful for this data. There is no competing alignment for source-quality to discount.
- The 71.6% G call at position 348 is **not** explained by paralog mapping. Whether it is a real germline het, an allelic-bias / PCR-amplification artifact, or sample-specific, is a biology question outside this document ;  but it can be looked at per-cell (see § 7).

## 5. Recommended option set (if LoFreq is used)

LoFreq has `-t/--approx-threshold INT` (default `-1`, off). Above a locus-depth threshold the exact DP is replaced by a Poisson-binomial to a Poisson approximation:

According to their source code:
```c
mu = sum err_probs[i];
pvalue = 1 - gsl_cdf_poisson_P(K - 1, mu);
```

Computation collapses to `O(N)` regardless of `K`. For Phred-quality error probabilities (~10⁻³), The approximation (Le Cam's theorem according to the source code) is accurate enough. The help-text in the code which says: "*might decrease number of calls*" applies only to borderline loci near the `sig / bonf` significance level, which do not exist at the depth and varcall boundary is `0.01 / 333 ≈ 3 × 10⁻⁵`; the calls in question have p-values at `10⁻⁴⁹³²`.

```
lofreq call \
  -f ref.fa \
  -o out.vcf.gz \
  --no-ext-baq \
  --no-default-filter \
  --min-bq 10 \
  --min-alt-bq 10 \
  --min-mq 50 \
  --max-mq 60 \
  --approx-threshold 10000 \
  aligned.sorted.bam
```

| Flag | Value | Rationale |
|---|---|---|
| `-e/--no-ext-baq` | on | Extended BAQ misbehaves at amplicon primer regions |
| `--no-default-filter` | on | Emit raw VCF; filter downstream for traceability |
| `-q/--min-bq` | `10` | Default `6` is too permissive |
| `-Q/--min-alt-bq` | `10` | Stricter on alt-supporting bases |
| `-m/--min-mq` | `50` | Justified by § 4.1 MAPQ histogram |
| `-M/--max-mq` | `60` | Explicit cap; prevents over-weighting outliers |
| `-t/--approx-threshold` | `10000` | Critical speed fix ;  `O(N × K)` → `O(N)` above this depth |

Deliberately not set, with reasons:

- **`-s/--src-qual`** ;  `s2:i:0` everywhere; no competing chains to discount (§ 4.3).
- **`-d/--max-depth`** ;  keep all signal; We should use `-t` not to lose coverage info.
- **`--call-indels`** ;  requires upstream `lofreq indelqual --dindel` to add IDAQ tags. Enable only if/when the indel preprocessing step is added.
- **`-C/--min-cov`** ;  irrelevant at this depth.
- **`--use-orphan`** ;  single-cell is R2 only barcode thus we don't need it.
- **`-a/--sig`**, **`-b/--bonf`** ;  defaults (`0.01`, `dynamic`) are correct.

For the Nextflow module the relevant entry is `task.ext.args` in [modules/lofreq/main.nf](../../modules/lofreq/main.nf).

## 6. Empirical verification of the options

Four LoFreq runs on the same BAM, compared with [`compare_caller_runs.sh`](compare_caller_runs.sh) (in this directory). Goal: confirm `-t` does what § 5 claims and characterise the parameter sensitivity around the recommended option set.

### 6.1 Runs

| Label | Flags (beyond `--bed`, `-f`, `--debug`, `--verbose`) | Threads | Wall | CPU (user) |
|---|---|---|---|---|
| `lofreq_100000t` | `-t 100000` (BAQ on, default filters) | 7 | 522m | 1,246m |
| `lofreq_100000t_250000d` | `-t 100000 -d 250000` (BAQ on, default filters) | 7 | 129m | 631m |
| `lofreq_10000t_minbq10_minmq50` | `-t 10000` + § 5 option set | 3 | 487m | 1,069m |
| `lofreq_100000t_minbq10_minmq50` | `-t 100000` + § 5 option set | 3 | 490m | 1,099m |

Thread counts vary across rows due to limitation running on AWS and limited time, so cross-row wall-time numbers mix algorithmic effects with parallelism unfortunately. The run time testing needs to be re-evaluated when/if we use LoFreq in production data.

### 6.2 Findings

**`-t 10000` vs `-t 100000` is not useful at this depth.** The two `_minbq10_minmq50` runs are identical: same 1,735 sites, Jaccard 1.0, max |dAF| = 0. At loci depth ~650k, both thresholds trigger the approximation on every loci. The § 5 recommendation of `-t 10000` stands and is valid.

**`--max-depth 250000` is the actual wall-time option once `-t` is engaged.** 522m to 129m (4×) at matched thread count. Computational Cost: 26 baseline sites drop and 8 new sites appear (597 shared, Jaccard 0.946), with one locus shifting VAF by ~0.22 (Pearson r = 0.998, median |dAF| = 0). Not added to § 5's recommended option set because it alter the call set, but it is the right next option if run time becomes binding.

**Quality filters dominate over `-t`.** Baseline (BAQ on, default filters) vs § 5 option set (`--no-ext-baq --min-bq 10 --min-alt-bq 10 --min-mq 50 --max-mq 60`): 623 sites vs 1,735 sites, 433 shared (Jaccard 0.225), median |dAF| 0.008 on the 413-site core. The two call sets are nearly un-identical, with § 5's option set recovering ~3× the calls ;  consistent with extended BAQ over-downgrading per-base qualities at amplicon primer regions (§ 5 rationale). LoFreq's call set is more sensitive to its filter configuration than to the algorithmic-approximation.

### 6.3 Scope

- Confirmed: `-t` option eliminates the § 3 stall, with bit-identical output across the `-t` operating range tested.
- Confirmed: § 5 option set is internally self-consistent.
- Open: LoFreq's call set against bcftools / freebayes / vardict / gatk-hc. The cross-caller comparison is what would either restore LoFreq or finalise exclusion. But the runtime of LoFreq does not make its inclusion attractive.

It is important to note that `compare_caller_runs.sh`'s produced Jaccard ≥ 0.99, max |dAF| < 0.02 are calibrated for same-parameter reproducibility, not cross-parameter comparison.

## 7. Decision

LoFreq is excluded from the default variant-caller set. The other implemented callers (bcftools, freebayes, vardict, gatk-hc) already cover germline SNV calling without the depth-dependent pathology. The `-t` mitigation is single-flag and accurate, but the behavior and the `-1` default make this failure mode tied to region read depth. This is a risk for even higher read depth data if that ever comes.

The implemented LoFreq module ([modules/lofreq/main.nf](../../modules/lofreq/main.nf)) is left in but it will be removed from the default varian callers in the subworkflow. Removing it is a separate issue.

## 8. Future notes

If a future developer wants to revive LoFreq or extend on this, here is what is open to do and discuss:

1. **Revive LoFreq, fast path:** add `task.ext.args = '--no-ext-baq --no-default-filter --min-bq 10 --min-alt-bq 10 --min-mq 50 --max-mq 60 --approx-threshold 10000'` to the LoFreq process in the config and add the LoFreq module into the default variant-calling subworkflow. § 6 confirms the option set runs without the depth stall (8.1h wall on the fixture) and is internally self-consistent across `-t`; the remaining gap is the cross-caller comparison against bcftools / freebayes / vardict / gatk-hc on the same fixture. If LoFreq does not produce calls the others miss, exclusion is even more fortified and easy decision.
2. **Per-cell VAF at position 348:** the 71.6% G call is a bulk sample. Stratify the same reads by cell barcode (`CB:Z:` or equivalent on the BAM) and recompute VAF per cell. Maybe a bimodal distribution (most cells het, some cells hom-alt) suggests a real germline variant plus maybe LOH or amplification effect; a continuous distribution suggests a PCR or amplification artifact. This is the biology question out of the question of this document; it is left for data analysis step.
3. **MAPQ-52 amplicon characterization:** § 4.3 explains the spike as a minimap2 `fac`-factor artifact but does not derive the exact reason. If a future minimap2 upgrade changes MAPQ distribution, the `--min-mq 50` threshold may need re-checking against a fresh histogram and some more data exploration (§ 4.1 command).
4. **Indels:** out of scope here. If indels are added to v1, `--call-indels` requires `lofreq indelqual --dindel` upstream of `lofreq call` there is no fix-up unfortunately. The calling step can do on its own.
5. **Source for the algorithm:** if the LoFreq DP needs re-reading, the central files are `src/lofreq/snpcaller.c` (`pruned_calc_prob_dist`, `poissbin`, `snpcaller`), `src/lofreq/lofreq_call.c` (`call_snvs`, where `alt_counts` / `noncons_counts` is populated), and `src/lofreq/defaults.h` (`NUM_NONCONS_BASES = 3`, `DEFAULT_SIG = 0.01`). The `-t` recommendation is in `snpcaller.c` is controlled by `conf->approx_threshold_n`.
