# NOTE THIS DOCUMENT IS HEAVILY EDITED BY AI AFTER DEVELOPMENT WAS DONE

# Variant Calling ; Technical Design

| Field | Value |
|---|---|
| Status | Draft |
| Scope | Germline SNV/SNP calling on amplicon BAMs |
| Out of scope (v1) | Indels, somatic mode, non-amplicon inputs |
| Touches | `main.nf`, `subworkflows/dna_processing.nf`, `modules/dna_processing.nf` (removals), new `subworkflows/variant_calling.nf`, per-tool `modules/<tool>/main.nf` (e.g. `modules/freebayes/main.nf`, `modules/samtools/faidx/main.nf`, `modules/bcftools/norm/main.nf`, `modules/bcftools/stats/main.nf`, `modules/tabix/main.nf`), `bin/`, `nextflow_schema.json`, `assets/schema_input.json`, `assets/metrics/metrics_reference.csv` |
| Migration | **Replacement, not additive.** The existing `PILEUP_READS` / `bcftools` pseudo-bulk path is removed; the new multi-caller subworkflow is the only variant-calling route. No feature flag. |

---

## Table of Contents

- [1. Goals and Non-Goals](#1-goals-and-non-goals)
  - [1.1 Goals](#11-goals)
  - [1.2 Non-Goals (v1)](#12-non-goals-v1)
- [2. Background](#2-background)
  - [2.1 Current variant-calling path (to be removed)](#21-current-variant-calling-path-to-be-removed)
  - [2.2 What's already in place (reuse targets)](#22-whats-already-in-place-reuse-targets)
  - [2.3 What's missing](#23-whats-missing)
- [3. Architecture](#3-architecture)
  - [3.1 Dataflow](#31-dataflow)
  - [3.2 Architecture and Channels in Nextflow](#32-architecture-and-channels-in-nextflow)
- [4. Module Specifications](#4-module-specifications)
  - [4.1 Process interface contract](#41-process-interface-contract)
  - [4.2 Caller CLIs (germline, SNV-only)](#42-caller-clis-germline-snv-only)
  - [4.3 Prepare reference inputs](#43-prepare-reference-inputs)
    - [4.3.1 `SAMTOOLS_FAIDX`](#431-samtools_faidx)
    - [4.3.2 `GATK4_CREATESEQUENCEDICTIONARY`](#432-gatk4_createsequencedictionary)
    - [4.3.3 `BEDTOOLS_MAKEWINDOWS`](#433-bedtools_makewindows)
  - [4.4 `BCFTOOLS_NORM`](#44-bcftools_norm)
  - [4.5 `BCFTOOLS_STATS`](#45-bcftools_stats)
  - [4.6 `VCF_TO_COUNTS_CSV` — REMOVED](#46-vcf_to_counts_csv--removed-collapsed-into-47)
  - [4.7 `VCF_TO_PARQUET` ; `bin/vcf_to_parquet.py`](#47-vcf_to_parquet--binvcf_to_parquetpy)
  - [4.8 `COUNT_VARIANT` (reused)](#48-count_variant-reused)
  - [4.9 `bin/variants_to_metrics_parquet.py`](#49-binvariants_to_metrics_parquetpy)
  - [4.10 `bin/qc_variants.py` (sibling to §4.9)](#410-binqc_variantspy-sibling-to-49)
- [5. Output Layout](#5-output-layout)
- [6. Configuration Reference](#6-configuration-reference)
  - [6.1 New pipeline params](#61-new-pipeline-params-nextflow_schemajson)
  - [6.2 Samplesheet additions](#62-samplesheet-additions-assetsschema_inputjson)
  - [6.3 Metrics reference](#63-metrics-reference-assetsmetricsmetrics_referencecsv)
- [7. Containers](#7-containers)
  - [7.1 Tool inventory](#71-tool-inventory)
  - [7.2 Creating and using containers](#72-creating-and-using-containers)
  - [7.3 Tool maintainability risk](#73-tool-maintainability-risk)
- [8. Testing Strategy](#8-testing-strategy)
  - [8.1 Unit / script-level](#81-unit--script-level)
  - [8.2 Integration (Nextflow)](#82-integration-nextflow)
  - [8.3 CI](#83-ci)
- [9. Implementation Plan and milestones](#9-implementation-plan-and-milestones)
- [10. Future Extensions](#10-future-extensions)
  - [10.1 Benchmarking](#101-benchmarking)
  - [10.2 Indel support](#102-indel-support)
  - [10.3 Somatic support](#103-somatic-support)
  - [10.4 Error correction](#104-error-correction)
- [11. Known Risks](#11-known-risks)
- [12. Open Decisions](#12-open-decisions)
- [Appendix A ; Files that will be affected (modify, create, delete)](#appendix-a--files-that-will-be-affected-modify-create-delete)
- [Appendix B ; Other existing workflows for variant callers](#appendix-b--other-existing-workflows-for-variant-callers)
- [Appendix C ; Cell barcode collapsing for error correction](#appendix-c--cell-barcode-collapsing-for-error-correction)

---

## 1. Goals and Non-Goals

### 1.1 Goals
1. Replace the "quite poor" pseudo-bulk `bcftools` SNP caller (see [README](README.md)) with a multi-caller bulk variant-calling path that runs up to **six callers in parallel** (freebayes, LoFreq, bcftools, VarDict, **DeepVariant**, **GATK HaplotypeCaller**) per sample. LoFreq is currently dormant in the default `variant_callers` — wired but excluded pending an INFO→FORMAT lift; the other five run by default.
2. Preserve the existing per-cell attribution logic in [bin/count_variant.py](bin/count_variant.py) by feeding it VCF-derived input that matches its current CSV schema.
3. Produce per-sample VCFs (bgzipped + tabix) and tidy long-format parquet tables with the following keys `(sample_name, caller)`.
4. Retire the legacy `PILEUP_READS` pseudo-bulk process entirely. This includes its downstream CSV cleanly. Any output previously consumed from the legacy path is either re-derived from the new VCFs or removed.

### 1.2 Non-Goals (v1)
- **Indel calling.** Callers are configured to skip indels; the normalization step filters to SNV records only. See §10 for the indel extension plan.
- **Somatic mode.** All callers run germline-configured. See §10 for the somatic extension plan.
- **Caller-level tumor/normal pairing.** v1 is single-sample.
- **Non-amplicon inputs** (WGS/WES).
- **Changes to `count_variant.py`.** The counter is SNV-only today and stays SNV-only in v1.

---

## 2. Background

uniflow is a Nextflow DSL2 pipeline (v0.4.14, alpha) that processes single-cell RNA+DNA co-sequencing libraries. Containers are Seqera Wave / Docker. The DNA modality is amplicon-targeted.

### 2.1 Current variant-calling path (**to be removed**)
`DNA_PROCESSING` ([subworkflows/dna_processing.nf](subworkflows/dna_processing.nf)) runs two branches off `ALIGN_AMPLICON.out.bam`:
- **Bulk SNP discovery:** `PILEUP_READS` ([modules/dna_processing.nf:77-125](modules/dna_processing.nf#L77-L125)) runs `bcftools mpileup --max-depth 10000 --skip-indels` followed by `call -mv` followed by `view --max-alleles 2` followed by CSV via `bcftools query`.
- **Per-cell attribution(s):** `COUNT_VARIANT` ([subworkflows/dna_processing.nf:42-49](subworkflows/dna_processing.nf#L42-L49)) joins the bulk CSV with the aligned BAM, barcode parquet, and amplicon metadata; `bin/count_variant.py` runs pysam pileups per cell. The legacy hardcoded `MIN_BASE_QUALITY=30` / `MIN_DP=200` thresholds were removed in COMB-436: site-level DP gating moved upstream to `BCFTOOLS_FILTER` (§4.4.1, driven by `params.variant_min_dp`), base-quality filtering moved to the variant callers themselves (§4.2 "Base-quality filtering"), and low-VAF noise is gated per-cell via `params.variant_min_vaf` inside `filter_snps()`. **Kept**, but rewired to consume VCF-derived input from the new subworkflow (see §4.8).

### 2.2 What's already in place (reuse targets)
- `bcftools 1.21` Wave image: `community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11`
- `samtools` / `tabix` bundled in minimap2+samtools and bcftools images
- Parallel BAM sharding utilities; debatable how this would fair when we start to scale for higher read counts: [bin/lib/bam_tools/utils.py](bin/lib/bam_tools/utils.py)
- Resource labels: `cpu_low_mem`, `small_job`, `lazy_polars`
- MultiQC mixing pattern: [main.nf:110](main.nf#L110) (FALCO); metrics YAML pattern at [subworkflows/dna_processing.nf:94](subworkflows/dna_processing.nf#L94)

### 2.3 What's missing
freebayes, LoFreq, VarDict, GATK ; not containerized in-repo. v1 sources these directly from BioContainers' pinned-hash images (`quay.io/biocontainers/<tool>:<version>--<build>`); see §7 for the policy and rationale. DeepVariant has an official Docker image from Google (`google/deepvariant`) and is consumed directly.

---

## 3. Architecture

### 3.1 Dataflow

Input: `ALIGN_AMPLICON.out.bam` → `VARIANT_CALLING` (sole variant-calling path; `PILEUP_READS` removed)

Variant calling will run per sample. A DNA sample runs through `VARIANT_CALLING` if:

```
run = sample.run_variant_calling ?? params.run_variant_calling
```

When `run` is false, the sample's BAM never enters `VARIANT_CALLING`; no caller, normalize, parquet, or `COUNT_VARIANT` step executes for it. The rest of `DNA_PROCESSING` (alignment, amplicon counting, QC) is unaffected. See §6.1, §6.2.

Variant callers (run in parallel, each will produce a `(meta, caller, vcf.gz)` tuple). Indexing (`.tbi`) and `bcftools stats` are factored out of caller modules and run as separate downstream steps (§4.1, §4.5).

BED dependency: `CALL_VARDICT` and `CALL_GATKHC` (M2), and `CALL_DEEPVARIANT` (M3), require an amplicon target BED. When `params.amplicon_bed` is null, the BED is derived once from the amplicon FASTA via `BEDTOOLS_MAKEWINDOWS` (§4.3.3) and broadcast to every BED-consuming caller; when `params.amplicon_bed` is set, the file is used directly and `BEDTOOLS_MAKEWINDOWS` does not run. `CALL_FREEBAYES`, `CALL_LOFREQ`, and `CALL_BCFTOOLS` do not consume the BED.

`CALL_FREEBAYES`
`CALL_LOFREQ`
`CALL_BCFTOOLS`
`CALL_VARDICT`
`CALL_DEEPVARIANT`
`CALL_GATKHC`
Normalization: All caller outputs flow into `BCFTOOLS_NORM` (bcftools norm; SNV filter; `QUAL/DP` gate)

Downstream outputs from `BCFTOOLS_NORM`:

`BCFTOOLS_STATS` → `MultiQC`
Per-caller VCFs → `VARIANT_HARMONIZE` → `VARIANT_MERGING` → merged consensus VCF

Downstream of the merged consensus VCF (one per sample):

`VCF_TO_PARQUET` → tidy long-format parquet, doubles as the per-cell catalog source
`PER_CELL_VARIANTS` → `COUNT_VARIANT` (per caller, filters the parquet by caller)

### 3.2 Architecture and Channels in Nextflow

A single entry point is designed to host all the downstream of the caller modules via `(meta, caller)` where `caller is a {freebayes, lofreq, bcftools, vardict, deepvariant, gatkhc}`. 

**IMPORTANT** Join across caller streams must group by this two-element `(meta, caller)` , never by `meta` alone.

The step for `DNA_PROCESSING to VARIANT_CALLING` condition, will not be inside the subworkflow. `DNA_PROCESSING` filters the BAM stream before invoking `VARIANT_CALLING`. A rough code sketch could be like:

```groovy
// subworkflows/dna_processing.nf
ALIGN_AMPLICON.out.bam
    | filter { meta, _bam, _bai ->
        meta.run_variant_calling != null ? meta.run_variant_calling : params.run_variant_calling
      }
    | set { vc_inputs }

VARIANT_CALLING(vc_inputs, amplicon_fasta, amplicon_fasta_fai, amplicon_fasta_dict, amplicon_bed)
```

This keeps the subworkflow itself unconditional and avoids leaking the workflow branching and logic into every downstream channels.

```groovy
// subworkflows/variant_calling.nf ; sketch
workflow VARIANT_CALLING {
    take:
        bam_with_index            // tuple(meta, bam, bai)
        amplicon_fasta            // path

    main:
        SAMTOOLS_FAIDX(amplicon_fasta)
        // GATK4_CREATESEQUENCEDICTIONARY lands in M2 alongside CALL_GATKHC
        GATK4_CREATESEQUENCEDICTIONARY(amplicon_fasta)

        // BED is conditional: external file wins, otherwise derive from the FASTA.
        // Kept inside VARIANT_CALLING so the subworkflow stays self-contained for
        // its reference inputs; DNA_PROCESSING does not need to know about BEDs.
        def amplicon_bed_ch
        if (params.amplicon_bed) {
            amplicon_bed_ch = channel.value(file(params.amplicon_bed, checkIfExists: true))
        } else {
            BEDTOOLS_MAKEWINDOWS(SAMTOOLS_FAIDX.out.fai)
            amplicon_bed_ch = BEDTOOLS_MAKEWINDOWS.out.bed
        }

        def caller_channels = []
        if ("freebayes"    in params.variant_callers) caller_channels << CALL_FREEBAYES(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai)
        if ("lofreq"       in params.variant_callers) caller_channels << CALL_LOFREQ(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai)
        if ("bcftools"     in params.variant_callers) caller_channels << CALL_BCFTOOLS(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai)
        if ("vardict"      in params.variant_callers) caller_channels << CALL_VARDICT(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai, amplicon_bed_ch)
        if ("deepvariant"  in params.variant_callers) caller_channels << CALL_DEEPVARIANT(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai, amplicon_bed_ch)
        if ("gatkhc"       in params.variant_callers) caller_channels << CALL_GATKHC(bam_with_index, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai, GATK4_CREATESEQUENCEDICTIONARY.out.dict, amplicon_bed_ch)

        raw_vcfs = Channel.empty().mix(*caller_channels.collect { it.vcf })
        BCFTOOLS_NORM(raw_vcfs, SAMTOOLS_FAIDX.out.fasta, SAMTOOLS_FAIDX.out.fai)

        // fan-outs from normalized stream
        TABIX(BCFTOOLS_NORM.out.vcf)
        BCFTOOLS_STATS(BCFTOOLS_NORM.out.vcf)

    emit:
        vcf   = BCFTOOLS_NORM.out.vcf
        tbi   = TABIX.out.tbi
        stats = BCFTOOLS_STATS.out.stats
}
```

`VCF_TO_PARQUET` and `COUNT_VARIANT` are not children of `BCFTOOLS_NORM` ; they consume the **merged consensus VCF** emitted from `VARIANT_MERGING.out.consensus_vcf` and live in [subworkflows/per_cell_variants.nf](subworkflows/per_cell_variants.nf), invoked from [subworkflows/dna.nf](subworkflows/dna.nf). See §4.7 and §4.8.

Why `BEDTOOLS_MAKEWINDOWS` lives inside `VARIANT_CALLING` rather than in `DNA_PROCESSING`: the BED is exclusively a variant-calling input (no other DNA step consumes it), and the conditional logic should live next to the consumers. `DNA_PROCESSING` continues to invoke `VARIANT_CALLING(vc_inputs, amplicon_fasta_file_ch)` with no BED-aware plumbing.

`DNA_PROCESSING` invokes `VARIANT_CALLING` as its sole bulk variant-calling branch. The `PILEUP_READS` process and all its content will be deleted Milestone 1 (See §9)

---

## 4. Module Specifications

### 4.1 Process interface contract

Since we are dealing more or less similar type of behaviour in variant callers, it is rationale to design modules in a way that every caller process has the same I/O shape:

```groovy
process CALL_<NAME> {
    label 'cpu_low_mem'
    container "..."  // see §7

    input:
        tuple val(meta), path(bam), path(bai)
        path fasta
        path fasta_fai
        // VarDict only:
        path amplicon_bed

    output:
        tuple val(meta), val("<caller>"), path("${meta.id}.<caller>.vcf.gz"), emit: vcf

    script:
    def args = task.ext.args ?: ""
    // TODO: design logging patterns and output path
    """
    <caller cmd> ${args} ...
    bgzip -f ${meta.id}.<caller>.vcf
    """
}
```

`task.ext.args` is pinned in `conf/modules.config`. v1 has no per-sample mode switching; every sample gets the same germline SNV args and the same treatment. `bgzip` stays in the caller module (it is the canonical compression step paired with the caller's text-VCF output), but `tabix` indexing and `bcftools stats` are factored out into standalone `TABIX` (`modules/tabix/main.nf`) and `BCFTOOLS_STATS` (§4.5) modules so caller modules stay tool-pure.

### 4.2 Caller CLIs (germline, SNV-only)

The following commands and args were taken from respective tools' repositories:

| Caller | Command sketch |
|---|---|
| freebayes | `freebayes --ploidy 2 --min-alternate-fraction 0.2 --min-alternate-count 2 --min-coverage 20 --no-indels --no-mnps --no-complex -f $fasta $bam` |
| LoFreq | `lofreq call-parallel --pp-threads $task.cpus -f $fasta -o out.vcf.gz $bam`  *(no `--call-indels`, no `indelqual` prereq; see "LoFreq notes" below)* |
| bcftools | `bcftools mpileup -a AD,DP,SP --skip-indels -f $fasta $bam \| bcftools call -mv --ploidy 2 -Oz -o out.vcf.gz` |
| VarDict | `VarDict -G $fasta -f 0.2 -r 2 -c 1 -S 2 -E 3 -g 4 -b $bam $amplicon_bed \| teststrandbias.R \| var2vcf_valid.pl -A -f 0.2 -N $sample > out.vcf` |
| DeepVariant | `run_deepvariant --model_type=WES --ref=$fasta --reads=$bam --regions=$amplicon_bed --output_vcf=out.vcf.gz --num_shards=$task.cpus` |
| GATK_HC | `gatk HaplotypeCaller -R $fasta -I $bam -L $amplicon_bed -O out.vcf.gz --output-mode EMIT_VARIANTS_ONLY --sample-ploidy 2` |

**Base-quality filtering**

Per-base quality filtering is delegated to the variant callers themselves; the per-cell pysam pileup in [bin/lib/modality/dna/processing.py](bin/lib/modality/dna/processing.py) (`pileup_truncated`) is an *attribution* step running with `min_base_quality=0`. Re-applying a base-quality cutoff at attribution time would discard reads the caller already accepted, so the responsibility lives upstream:

- freebayes — `--min-base-quality` (default 0; raise with caller args if needed)
- bcftools — `bcftools mpileup -Q` (default 13)
- LoFreq — `--min-bq` / `--min-alt-bq` (default 6 / 6)
- VarDict — `-q` (default 22.5)
- DeepVariant — handled internally by the model (no exposed knob)
- GATK HaplotypeCaller — `--base-quality-score-threshold` (default 18)

Tuning per-caller base-quality thresholds is out of scope for v1; defaults stand until a truth-set comparison drives a change.

**LoFreq notes**
- LoFreq's VCF output **omits the `##contig=<ID=,length=>` header lines** that other callers propagate from the BAM's `@SQ` records. Downstream consumers (notably GATK-family tools, but also `bcftools view -h` consistency checks) expect them. To fix this in-place without touching the LOFREQ module, the lofreq stream goes through `GATK4_UPDATEVCFSEQUENCEDICTIONARY` (`gatk UpdateVCFSequenceDictionary --source-dictionary <bam> --replace true --create-output-variant-index false`) before flowing into `BCFTOOLS_NORM`. The BAM is used as the dictionary source rather than the `.dict` from `GATK4_CREATESEQUENCEDICTIONARY` because (a) the BAM is the authoritative record of what contigs lofreq actually called against, and (b) it removes the dependency on the `.dict` module, which is only built when `gatkhc` is in `params.variant_callers`.
- `UpdateVCFSequenceDictionary` requires an indexed input VCF; the module runs `gatk IndexFeatureFile --input <vcf>` as its first step. Indexing inside the reheader module (rather than emitting a `.tbi` channel from LOFREQ) keeps LOFREQ's output shape identical to the other callers and confines the lofreq-specific quirk to one place. The cost is one extra JVM cold-start per lofreq sample.
- The reheadered VCF is what enters `caller_vcfs`, so from `BCFTOOLS_NORM` onward lofreq is handled identically to every other caller.

**VarDict notes**
- Implemented in [modules/vardict/main.nf](modules/vardict/main.nf) (COMB-418). Uses `vardict-java` (1.8.3). It comes with `teststrandbias.R` and `var2vcf_valid.pl`.
- `var2vcf_valid.pl` is single-sample (germline). `var2vcf_paired.pl` is reserved for the somatic extension (§10).
- The BioContainers `vardict-java:1.8.3--hdfd78af_0` image does not bundle `bgzip`, so the module emits an uncompressed `.vcf` and the subworkflow runs a separate `BGZIP` step ([modules/tabix/bgzip/main.nf](modules/tabix/bgzip/main.nf), htslib via the bcftools image) before joining `caller_vcfs`. Output of the bgzip step is `${meta.sample_name}.vardict.vcf.gz`, matching the channel shape of the other callers.
- Target BED is mandatory. Fallback: derive from amplicon FASTA by treating each record as a single interval (`chrom 0 len`). Put this derivation in `BEDTOOLS_MAKEWINDOWS` (enable/disable it on `params.amplicon_bed == null`). Implemented via `bedtools makewindows -g <fai> -n 1` ; see §4.3.3.
- We drop any indel records in `BCFTOOLS_NORM` regardless of `-I 0` behavior. This behavior will change in expansion phase, see §10.

**DeepVariant notes**
- The official `google/deepvariant` Docker image (no Wave rebuild ; see §7).
- `--model_type=WES` is the closest match for short-read amplicon panels. WGS is a reasonable alternative; the WES model has been trained on more targeted-enrichment-like coverage distributions, which is closer to amplicon than WGS. Revisit if performance underperforms on the amplicon fixture.
- `--regions=$amplicon_bed` is non-optional here ; without it, DeepVariant examines the entire reference and wastes most of its runtime on zero-coverage regions of the amplicon FASTA. Uses the same `BEDTOOLS_MAKEWINDOWS` fallback as VarDict.
- DeepVariant will produce SNVs and indels natively; indels are stripped by `BCFTOOLS_NORM` in v1 (same as VarDict).
- Resource profile is heavier than the other callers (CPU-bound by the `make_examples` stage). `TODO: Consider "label 'cpu_mid_mem'" or a dedicated label.` GPU acceleration is out of scope for v1 and out of scope of the whole project for beta launch in September.

**GATK HaplotypeCaller notes**
- `gatk4` Docker image (see §7). Requires both a `.fai` and a `.dict` sidecar; the `.dict` is produced by `GATK4_CREATESEQUENCEDICTIONARY` (§4.3.2).
- `-L $amplicon_bed` restricts to the target panel. Same BED-fallback rule as VarDict/DeepVariant.
- `--output-mode EMIT_VARIANTS_ONLY` keeps the VCF to variant sites only (no GVCF). GVCF + joint genotyping is deferred; for single-sample amplicon the direct VCF mode is correct.
- HaplotypeCaller just like DeepVariant will produce SNVs and indels; indels are filtered in `BCFTOOLS_NORM` for v1.
- Java heap: set via `--java-options "-Xmx8g"` **hardcoded inside the module body** (not derived from `task.memory`). GATK's container defaults under-allocate, and the amplicon-scale workload does not need (or benefit from) the full `cpu_low_mem` 20 GB; pinning `-Xmx8g` keeps behaviour predictable across resource labels and leaves clear native/JVM-overhead headroom on `cpu_low_mem`. The constant lives in the module, not in `conf/modules.config`, so memory tuning is a one-line module edit. Revisit if amplicon panels grow large enough to spill the heap, or migrate to a dedicated `cpu_mid_mem` label at that point.

### 4.3 Prepare reference inputs

### 4.3.1 `SAMTOOLS_FAIDX`

A simple `samtools faidx` wrapper. Creates `.fai` next to the fasta. Pure convenience, butterflies and puppies; every caller needs it anyway. Lives in `modules/samtools/faidx/main.nf` so future `samtools` subcommands (e.g. `dict`, see §4.3.2) can land beside it under `modules/samtools/<subcommand>/`.

### 4.3.2 `GATK4_CREATESEQUENCEDICTIONARY`

A `gatk CreateSequenceDictionary` wrapper producing `${fasta.baseName}.dict`. Required by GATK HaplotypeCaller. Runs at most once per pipeline run and is stored alongside the `.fai`. Lives at `modules/gatk4/createsequencedictionary/main.nf` and reuses the GATK 4.6.1.0 BioContainers image already pinned for `CALL_GATKHC` (§4.2 / §7).

### 4.3.3 `BEDTOOLS_MAKEWINDOWS`

Derives a target BED from the amplicon FASTA when `params.amplicon_bed` is `null`. One BED record per FASTA contig, with `start = 0` and `end = <contig length>` ; treats each amplicon record as a single full-length interval. Required by `CALL_VARDICT` (M2), `CALL_GATKHC` (M2), and later `CALL_DEEPVARIANT` (M3). When `params.amplicon_bed` is set, the user-supplied BED is used directly and this module does not run.

- Lives at [modules/bedtools/makewindows/main.nf](modules/bedtools/makewindows/main.nf). Uses `bedtools makewindows` from the BioContainers image `quay.io/biocontainers/bedtools:2.31.1--h13024bc_3`. **New container** ; deviates from the original "no new container" plan because (a) `bedtools makewindows -n 1` is the semantic primitive for "one full-length interval per contig" and is more legible than an awk one-liner, and (b) it sets up cleanly for the indel-aware extension (§10.2) where per-contig padding/slop becomes desirable. The trade is one extra image pull per pipeline run when `params.amplicon_bed` is null.
- Input: the `${fasta}.fai` already produced by `SAMTOOLS_FAIDX` (§4.3.1). `bedtools makewindows -g` accepts the `.fai` directly because it only consumes the first two columns (chrom, length) ; this avoids a redundant `samtools faidx` invocation and keeps a single canonical contig-length source.
- Output: `${fasta.simpleName}.bed`, emitted as a value channel and broadcast to every BED-consuming caller.
- Implementation sketch:

  ```bash
  bedtools makewindows -g ${fasta}.fai -n 1 > ${fasta.simpleName}.bed
  ```

- Edge cases:
  - **Empty FASTA / `.fai`** → empty BED. VarDict/GATK_HC will produce zero-row VCFs; `BCFTOOLS_NORM` and downstream steps run cleanly on empty inputs. Acceptable; we do not error here.
  - **Contigs with non-ASCII names or whitespace** → `samtools faidx` already rejects these, so we do not re-validate.
  - **User-supplied `params.amplicon_bed` whose contig names disagree with the FASTA** → out of scope for this module; flagged at caller runtime via the tool's own error.

- Not parameterised. There is no per-amplicon padding, no slop, no contig allow-list ; v1 derives whole-record intervals only. If an indel-aware extension (§10.2) needs slop, add it as a separate module rather than overloading this one.

- The conditional sits inside `VARIANT_CALLING` (see §3.2 sketch), not in `DNA_PROCESSING`, because the BED is exclusively a variant-calling input.

### 4.4 `BCFTOOLS_NORM`

```bash
bcftools norm -f $fasta -m -both --rm-dup all $in.vcf.gz \
  | bcftools view -v snps -e 'ALT="*"' \
  -Oz -o $out.vcf.gz
```

- `bcftools norm -m -both` splits multi-allelics and left-aligns. Left-alignment will not affect SNVs but is kept for forward compatibility with the indel phase.
- `bcftools view -v snps` is the hard SNV-only filter. Any indel/MNV that leaks past caller flags is dropped here. See expansion §10 why this applies.
- **QUAL gating is intentionally not applied here.** Per-caller QUAL scales differ enough (LoFreq ranges in the hundreds, freebayes/bcftools in the tens, DeepVariant phred-probability) that a single cutoff is not defensible and per-caller cutoffs would need to be re-tuned against truth data we do not yet have. Variant confidence is instead surfaced downstream via `INFO/NCALLERS` from the consensus VCF and the per-cell VAF floor in `count_variant.py`. See [clonal_assignment_mitigations.md §5.2](clonal_assignment_mitigations.md#52-tier-0-proxy-from-already-available-signals) for the tier-0 confidence proxy this leans on.
- Site-level DP filtering is **not** in this module either; it lives in §4.4.1 `BCFTOOLS_FILTER` immediately downstream.
- Indexing (`.tbi`) is a separate `TABIX` step (`modules/tabix/main.nf`) chained after `BCFTOOLS_FILTER` in the subworkflow.

### 4.4.1 `BCFTOOLS_FILTER`

```bash
bcftools view -e 'INFO/DP<${params.variant_min_dp}' -Oz -o $out.vcf.gz $in.vcf.gz
```

- Site-level hard depth filter on the per-caller normalised VCF; records with `INFO/DP < params.variant_min_dp` are dropped before merge. See §6.1 for the param.
- Lives at [modules/bcftools/filter/main.nf](modules/bcftools/filter/main.nf). The filter expression is inlined in the script body (not configured via `conf/modules.config`/`ext.args`) so the depth gate is reviewable as part of the module diff and the module stays a one-purpose hard filter rather than a configurable `bcftools view` shim.
- All downstream consumers (`TABIX`, `BCFTOOLS_STATS`, `VARIANT_HARMONIZE` → `VARIANT_MERGING`, and `PER_CELL_VARIANTS`) read from `BCFTOOLS_FILTER.out.vcf`, not `BCFTOOLS_NORM.out.vcf`.

### 4.5 `BCFTOOLS_STATS`

`bcftools stats` with caller-tagged output filename so MultiQC can separate sections. This will make it easier for downstream aggregiation and reporting. Lives at `modules/bcftools/stats/main.nf`; consumes `(meta, caller, vcf)` from `BCFTOOLS_NORM`. Caller modules do not run `bcftools stats` themselves (§4.1).

**MultiQC wiring + publish-path move (first M2 ticket).** The stats files do two jobs simultaneously: feed MultiQC and surface as published per-sample metric artefacts. Both are wired off the same channel.

1. `VARIANT_CALLING` already emits `stats = BCFTOOLS_STATS.out.stats` shaped as `(meta, caller, path)`.
2. `DNA_PROCESSING` adds a dedicated emit ; `variant_stats = VARIANT_CALLING.out.stats` ; rather than re-using the publish-only `variant_outputs` group, so MultiQC sees only stats files (not VCFs/TBIs) and the publish path can decouple from the VCF layout.
3. `main.nf` does two things with the same channel:
   - **MultiQC mix**: `ch_multiqc_files = ch_multiqc_files.mix(DNA_PROCESSING.out.variant_stats.map { _meta, _caller, file -> file })` ; same pattern as FALCO ([main.nf:113](main.nf#L113)). Result is one MultiQC `bcftools stats` section per `(sample, caller)` ; five sections per sample after M2.
   - **Publish**: stats land at `samples/{sample}/metrics/variants/` (flat, no per-caller subdirectory). Filename `<sample>.<caller>.bcftools_stats.txt` already disambiguates callers and matches MultiQC's `bcftools` module auto-detection pattern, so no custom MultiQC config is needed and per-caller nesting would be redundant. The path is set in `main.nf`'s `output {}` block via `path { meta, caller, file -> "samples/${meta.sample_name}/metrics/variants/" }`.
4. **Migration.** This moves the freebayes stats file from `samples/{sample}/DNA/variants/freebayes/` (M1's location) to `samples/{sample}/metrics/variants/`. M1 is alpha and has no external consumers of this path, so the move is a one-shot break ; no compat shim. The same publish path applies forward to all M2/M3 callers.

### 4.6 `VCF_TO_COUNTS_CSV` — REMOVED (collapsed into §4.7)

The originally planned counts-CSV intermediate was dropped during M4 implementation. The long-format parquet from §4.7 is a strict superset of what the CSV would have carried (same rows, more columns, per-caller QUAL via `INFO/VARQUAL[idx]`), so a separate CSV step would only duplicate information. [bin/count_variant.py](bin/count_variant.py) reads the parquet directly via `pl.scan_parquet(...).filter(pl.col("caller") == caller_name)` ; no CSV ever exists on disk.

Tickets that referenced this module (COMB-432 spec, COMB-437 schema enhancement) close out as superseded by §4.7.

### 4.7 `VCF_TO_PARQUET` ; `bin/vcf_to_parquet.py`

Single-source-of-truth catalog: reads the **merged consensus VCF** ([subworkflows/variant_merging.nf](subworkflows/variant_merging.nf) `consensus_vcf` emit, one per sample, sample columns named after callers) and writes a long-format parquet.

| Column | Type | Source |
|---|---|---|
| `sample_name` | str | process val (meta.sample_name) |
| `caller` | str | sample column name in the merged VCF |
| `chrom` | str | `#CHROM` |
| `pos` | int64 | `POS` (1-based) |
| `ref` | str | `REF` |
| `alt` | str | `ALT` (single allele; `BCFTOOLS_NORM -m -both` split multi-allelics upstream) |
| `qual` | float64 | `INFO/VARQUAL[idx]` ; per-caller site QUAL preserved by `BCFTOOLS_FILL_TAGS`, joined in sample-column order by `BCFTOOLS_MERGE --info-rules VARQUAL:join`, kept through `BCFTOOLS_ANNOTATE_CONSENSUS` |
<<<<<<< HEAD
=======
| `filter` | str | `FILTER` (joined `;` if multiple) |
>>>>>>> remotes/origin/main
| `dp` | int64 | per-caller `FORMAT/DP` |
| `ad_ref` | int64 | per-caller `FORMAT/AD[0]` |
| `ad_alt` | int64 | per-caller `FORMAT/AD[1]` |
| `af` | float64 | per-caller `FORMAT/VAF` (computed post-merge from AD by `BCFTOOLS_ANNOTATE_CONSENSUS`) |
| `caller_specific` | str (JSON) | empty string `"{}"` in v1; deferred (no per-caller INFO is currently passed through the merge) |

<<<<<<< HEAD
One row per (sample, caller, variant) where the caller has `FORMAT/DP > 0` for that record (i.e. that caller contributed the call). Rows where a caller did not call the variant are skipped, not emitted with `dp=null`. Empty VCF → zero-row parquet with the full 12-column schema, never an error.
=======
One row per (sample, caller, variant) where the caller has `FORMAT/DP > 0` for that record (i.e. that caller contributed the call). Rows where a caller did not call the variant are skipped, not emitted with `dp=null`. Empty VCF → zero-row parquet with the full 13-column schema, never an error.
>>>>>>> remotes/origin/main

Module: [modules/vcf_to_parquet/main.nf](modules/vcf_to_parquet/main.nf). Input shape: `(meta, vcf)`. Output: `(meta, parquet)`. Container reuses `community.wave.seqera.io/library/pip_pyyaml_upsetplot_numpy_pruned:3f7a3b0ec65284cc` (polars + pysam). One process invocation per sample, not per caller — the per-caller dimension is encoded as the `caller` column.

### 4.8 `COUNT_VARIANT` (reused, per caller)

Runs once per `(meta, caller)`. Input tuple: `(meta, caller, bam, barcode_parquet, fasta, calls_parquet, bai, filtered_amplicon_reads, anchors_parquet)`. The `calls_parquet` is the **same** per-sample parquet emitted by §4.7, staged into every per-caller task; [bin/count_variant.py](bin/count_variant.py) filters it on `caller == sys.argv[8]` and projects to `(CHROM, POS, REF, ALT, DP, QUAL)` before the existing pileup/attribution pipeline. Fan-out lives in [subworkflows/per_cell_variants.nf](subworkflows/per_cell_variants.nf) via `VCF_TO_PARQUET.out.parquet.combine(channel.fromList(callers))` + sequential `combine(..., by: 0)` for the static per-sample inputs. Publishing is namespaced by caller (see §5).

### 4.9 `bin/variants_to_metrics_parquet.py`

Reads the per-sample long-format parquet from `VCF_TO_PARQUET` (§4.7) and emits the metrics YAML that `DNA_PROCESSING.out.metrics_parquet` consumes. Mixing pattern matches [`COUNT_AMPLICON.out.metrics_parquet`](subworkflows/dna_processing.nf#L94).

The parquet alone is sufficient — by §4.7 construction it has one row per `(sample, caller, variant)` where the caller has `FORMAT/DP > 0`, so grouping by `(chrom, pos, ref, alt)` and counting distinct contributing callers is equivalent to `INFO/NCALLERS` from the merged consensus VCF. Reading the consensus VCF in addition would only re-derive the same signal and risk parquet/VCF drift.

Active caller set as of M5 design: `{freebayes, bcftools, vardict, gatkhc}` (lofreq dormant pending INFO→FORMAT lift, deepvariant lands in M3). YAML output scales by which callers actually contributed in the parquet — callers with zero contributions emit no rows (callers in `params.variant_callers` that never produced a call don't get zero-filled metrics in v1; if this causes cross-sample aggregation gaps it can be revisited).

**YAML output (flat keys).** `MetricsRecord` ([bin/lib/metrics/metrics.py](bin/lib/metrics/metrics.py)) stores metrics as a flat `Dict[str, Any]`, and `MetricsCollector.generate_df` skips non-numeric values — nested-dict structures would be silently dropped on aggregation. Keys are therefore `{caller}.{metric}` / `consensus.{metric}`:

```yaml
source_id: S1
library_type: DNA
modality: variants
process_name: variants_to_metrics_parquet
metrics:
  freebayes.snv_count: 1234
<<<<<<< HEAD
=======
  freebayes.pass_count: 1180
>>>>>>> remotes/origin/main
  freebayes.singleton_count: 87
  freebayes.mean_dp: 412.3
  freebayes.median_dp: 380.0
  freebayes.mean_vaf: 0.34
  freebayes.median_vaf: 0.28
  freebayes.het_count: 920
  freebayes.hom_count: 260
  bcftools.snv_count: ...
  vardict.snv_count: ...
  gatkhc.snv_count: ...
  consensus.total_sites: 1307
  consensus.high_confidence_count: 1212
  consensus.singleton_count: 95
  consensus.fraction_high_confidence: 0.927
```

**Metric definitions** (descriptions match the `description` column of [assets/metrics/metrics_reference.csv](assets/metrics/metrics_reference.csv) added in COMB-467):

| Metric | Source | Description |
|---|---|---|
| `<caller>.snv_count` | parquet rows filtered to `caller==X` | Normalized SNV sites this caller contributed |
<<<<<<< HEAD
=======
| `<caller>.pass_count` | `filter == "PASS"` rows | Caller-PASS sites. Caller-specific FILTER semantics — not cross-comparable |
>>>>>>> remotes/origin/main
| `<caller>.singleton_count` | parquet sites with `ncallers==1` AND this caller contributed | Sites this caller called alone |
| `<caller>.mean_dp` / `median_dp` | per-caller `FORMAT/DP` column | Per-site read depth |
| `<caller>.mean_vaf` / `median_vaf` | per-caller `FORMAT/VAF` column | Per-site variant allele fraction |
| `<caller>.het_count` | per-caller `FORMAT/VAF` in `(0.05, 0.8)` | Heterozygous-leaning sites |
| `<caller>.hom_count` | per-caller `FORMAT/VAF` ≥ 0.8 | Homozygous-leaning sites |
| `consensus.total_sites` | distinct `(chrom, pos, ref, alt)` in parquet | Total normalized SNV sites in the merged consensus VCF |
| `consensus.high_confidence_count` | sites with `ncallers >= 2` | Sites called by ≥2 callers. **Absolute threshold, not relative to the active caller count** — running with only 2 callers turns this into "called by all", a different bar. |
| `consensus.singleton_count` | sites with `ncallers == 1` | Inverse of high-confidence |
| `consensus.fraction_high_confidence` | `high_confidence_count / total_sites` | Robust to total-variant-count drift; preferred over the absolute count for cross-sample comparison |

**Deferred to a follow-up (not in this PR):**
- `per_cell.<caller>.{cells_with_calls,mean_variants_per_cell,median_variants_per_cell}` — needs `COUNT_VARIANT`'s per-caller h5ads (one per `(meta, caller)`) joined back to a sample-level group before invoking the metrics script. Skipped from v1 to keep the M5 first cut tight; revisit once the per-caller surface looks right in MultiQC.
- `consensus.ncallers_histogram` — dict-valued, dropped by `MetricsCollector.generate_df`. Could be flattened to `consensus.ncallers_1..N` but the UpSet plot (§4.10) already carries the agreement-distribution signal visually; redundant for v1.

**`report` column mapping (COMB-467).** Per-caller and `consensus.*` metrics → `aggregated` (cross-sample comparability matters).

**Thresholds (COMB-467).** All `min_warning_threshold` / `max_warning_threshold` / `min_failure_threshold` / `max_failure_threshold` columns left **blank** in v1. Matches the existing convention for uncalibrated metrics (`total_reads`, `filtered_cells`, etc.). Populate once benchmarking (§10.1) lands or after the first production rollout calibrates expected ranges.

**Process shape.** Module `VARIANTS_TO_METRICS_YAML` ([modules/variants_to_metrics_parquet/main.nf](modules/variants_to_metrics_parquet/main.nf)) takes `(meta, parquet)` and emits `(meta, "metrics.parquet")` — same emit shape as `COUNT_AMPLICON.out.metrics_parquet` ([modules/dna_processing.nf:19](modules/dna_processing.nf#L19)). Mixed into the existing `metrics_parquet` channel at [subworkflows/dna.nf:101](subworkflows/dna.nf#L101). Reuses the `pip_pyyaml_upsetplot_numpy_pruned` container (polars already present; pysam not needed since the script reads the parquet directly).

### 4.10 `bin/qc_variants.py` (sibling to §4.9)

Visual companion to the YAML metrics. Sibling script in the spirit of `qc_amplicon.py` paired with `count_amplicon.py`: separates "metric extraction" from "QC plotting" so the two evolve independently.

**Primary output: caller-agreement UpSet plot.**
- **Sets** = active callers in `params.variant_callers`
- **Elements** = normalized SNV sites, keyed `(chrom, pos, ref, alt)`
- **Membership** = caller contributed a call at this site (any row in the per-sample parquet for that `(caller, site)`)
- **Bars** = count of variants per intersection (called-by-`gatkhc`-only, called-by-all-four, called-by-`{freebayes,vardict}`, etc.)

Reuses `lib.plotting.common.plot_generic_upset_plot` ([bin/qc_amplicon.py:5](bin/qc_amplicon.py#L5)) and the existing `pip_pyyaml_upsetplot_numpy_pruned` container — no new container, no new dependency.

**Plot decisions (locked):**
- **Highlight: `NCALLERS >= 2` in green.** The high-confidence intersections (the same cutoff used by `consensus.high_confidence_count` in §4.9) are the visually-emphasised tier. Analogous to how `qc_amplicon.py` highlights the "Usable" intersection.
- **Data scope: all rows in the per-sample parquet.** No `FILTER` filter, no VAF filter, no `min_subset_size` cutoff. The parquet has already been through `BCFTOOLS_NORM` (SNV-only) and `BCFTOOLS_FILTER` (`INFO/DP >= params.variant_min_dp`, default 50); past that point we want the unfiltered view in v1.
- **Per-caller `FILTER` (`PASS` / non-`PASS`) is intentionally not stratified.** Each caller's `PASS` semantics differ enough that mixing them in one plot would compare apples to oranges.
- **`min_pct_show = None` (no minimum subset size).** Plot every intersection; revisit if it turns out noisy in real-data runs.

**What the UpSet does and does not see.**
- `params.variant_min_vaf` (default `0.01`) is applied **per-cell** inside [bin/count_variant.py](bin/count_variant.py) (§6.1), downstream of the parquet. It does not affect the parquet rows the UpSet reads. A future reader expecting `variant_min_vaf` to influence the plot would be surprised — it doesn't.
- `params.variant_min_dp` **does** affect the parquet (and therefore the UpSet) via `BCFTOOLS_FILTER` (§4.4.1). Sites below the depth floor are dropped before they reach the parquet.

**Secondary output: intersection-counts YAML.** Same data behind the UpSet bars, persisted as a flat `MetricsRecord` for cross-sample aggregation. Keys take the form `upset.<sorted_caller_intersection>` (e.g. `upset.freebayes_gatkhc.sites` for "called only by freebayes and gatkhc"). Useful for "did the caller-agreement pattern shift between batches?" without re-opening the PNG. Emitted alongside the PNG; mixed into `DNA.out.metrics_parquet` the same way `VARIANTS_TO_METRICS_YAML.out.metrics_parquet` is.

**Process shape.** New `QC_VARIANTS` module at [modules/qc/qc_variants/main.nf](modules/qc/qc_variants/main.nf), sibling to `modules/qc/caller_compare/main.nf` from M4. Takes `(meta, parquet)`, emits:
- `png = (meta, path("caller_upset.png"))` — canonical visual artifact, published to `samples/{sample}/qc/variants/caller_upset.png`, alongside the `pairwise_jaccard.tsv` / `summary.txt` already produced there by `CALLER_COMPARE_QC`.
- `png_mqc = (meta, path("caller_upset_mqc.png"))` — a `cp` of the canonical PNG renamed to satisfy MultiQC's `*_mqc.{png,jpg}` auto-discovery. Mixed into `ch_multiqc_files` so the UpSet appears as a per-sample image section in the aggregated report.
- `metrics_parquet = (meta.sample_name, path("metrics.parquet"))` — intersection counts as a `MetricsRecord` with flat `upset.<sorted_callers>` keys.

**Why two PNG outputs:** MultiQC's image-section auto-discovery requires a `*_mqc.*` filename suffix. Renaming the canonical artifact to `caller_upset_mqc.png` everywhere would leak a MultiQC-specific naming convention into the published output tree. A `cp` produces both files with one plot generation pass — cheaper than templating a YAML sidecar (`*_mqc.yaml` with `image:` content) and keeps the canonical name clean.

---

## 5. Output Layout

The following output structure is considering the whole implementation Milestones (specifically Milestone 4):

```
samples/{sample}/DNA/variants/{caller}/
    {sample}.{caller}.vcf.gz
    {sample}.{caller}.vcf.gz.tbi
    {sample}.{caller}.variants.parquet
    # if we have params.raw_vcf = true:
    {sample}.{caller}.raw.vcf.gz

samples/{sample}/DNA/variants/{caller}/per_cell/
    variants_count.h5ad
    variants_count.parquet
    raw_variants.parquet

samples/{sample}/metrics/variants/
    {sample}.{caller}.bcftools_stats.txt    # one per (sample, caller); flat, no per-caller subdirectory
    metrics.parquet                         # M5: §4.9 per-caller + consensus + per-cell

samples/{sample}/qc/variants/
    caller_upset.png                        # M5: §4.10 caller-agreement UpSet
    # plus the CALLER_COMPARE_QC artifacts from M4:
    pairwise_jaccard.tsv
    af_pairwise_stats.tsv
    summary.txt
```

`bcftools stats` files live under `metrics/variants/` rather than `DNA/variants/{caller}/`. Rationale: they are MultiQC-consumed metric artefacts, not variant-data outputs ; co-locating with the other per-sample metrics (`metrics/amplicon`, `metrics/gene_expression`, `metrics/carryover`) keeps the metric surface uniform. The flat layout is intentional ; `<sample>.<caller>.bcftools_stats.txt` already disambiguates callers, MultiQC auto-discovers them by filename, and per-caller subdirectories would just add a layer with no information.

`raw_vcf` (boolean, default `true`) controls whether the pre-normalization caller output is generated for debugging the SNV filter and for the future indel phase.

---

## 6. Configuration Reference

### 6.1 New pipeline params ([nextflow_schema.json](nextflow_schema.json))

| Param | Type | Default | Description |
|---|---|---|---|
| `run_variant_calling` | boolean | `true` | Pipeline-wide condition. When `false`, no DNA sample runs `VARIANT_CALLING`. Per-sample samplesheet column overrides this (§6.2). |
| `variant_callers` | array<string> | `["freebayes", "lofreq", "bcftools", "vardict", "deepvariant", "gatkhc"]` | Which callers to run. Subset is permitted. `minItems: 1` ; an empty array will raise a schema error and instruct the user in order to disable variant calling, use `run_variant_calling: false` instead. |
| `variant_min_dp` | integer | `50` | Site-level `INFO/DP` floor applied in `BCFTOOLS_FILTER` (§4.4.1) downstream of `BCFTOOLS_NORM`. Records below the floor are dropped from the per-caller VCFs before merge. Replaces the legacy hardcoded `MIN_DP=200` in `count_variant.py`; per-cell depth gating is no longer applied there (the parquet's per-caller `FORMAT/DP` is a caller-side report, not a per-cell count). |
| `variant_min_vaf` | number | `0.01` | Per-cell hard floor on per-variant VAF (ALT reads / total reads at that variant in that cell), applied in [bin/count_variant.py](bin/count_variant.py) inside `filter_snps()` before the homo/hete classification. Intended to suppress low-VAF sequencing/PCR/mapping errors that no upstream caller filter can reach. Does **not** affect the homo/hete labels of rows that pass the floor — `FRACTION_READS_HOMO=0.8` is unchanged. |
| `amplicon_bed` | file-path | `null` | Target BED for VarDict, DeepVariant, and GATK_HC. When `null`, derived once from `amplicon_fasta` via `BEDTOOLS_MAKEWINDOWS` (§4.3.3) ; one record per FASTA contig, full-length interval. When set, the user-supplied BED is used directly and the fallback module does not run. |
| `deepvariant_model_type` | string enum | `"WES"` | One of `WGS` \| `WES` \| `PACBIO` \| `ONT_R104` \| `HYBRID_PACBIO_ILLUMINA`. Amplicon panels use `WES`. |
| `raw_vcf` | boolean | `true` | Publish pre-normalization caller output alongside the normalized VCF. |

Note: `variant_mode` is not yet considered here. It is reserved for the somatic extension (see §10).

### 6.2 Samplesheet additions ([assets/schema_input.json](assets/schema_input.json))

| Column | Type | Required | Notes |
|---|---|---|---|
| `run_variant_calling` | boolean | No | Per-sample override of `params.run_variant_calling`. Empty cell means "use the global default". Only meaningful for `library_type == DNA` rows; ignored otherwise. |

### 6.3 Metrics reference ([assets/metrics/metrics_reference.csv](assets/metrics/metrics_reference.csv))

Add one row per `(caller, metric)` combination. Follow the `variants.{caller}.{metric}` convention (e.g., `variants.freebayes.snv_count`).

---

## 7. Containers

uniflow's existing container model is **Seqera Wave with pinned-hash URIs** for the pipeline's pre-variant-calling modules (e.g. `community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11`). Wave is *not* invoked at runtime; there is no `wave.enabled = true` in [nextflow.config](nextflow.config), and no module declares `conda::` directives.

**For the v1 variant-calling additions, we source BioContainers pinned-hash images directly** (`quay.io/biocontainers/<tool>:<version>--<build>`). Each tool listed in §7.1 below is a content-addressed BioContainers URI built upstream from the corresponding Bioconda recipe; we paste the URI into the module's `container "..."` line. This avoids the Wave-build round-trip, decouples our release cadence from Bioconda PR review of new Wave builds, and keeps reproducibility properties equivalent to the existing Wave-pinned model (both are content-addressed). See §7.2 for the policy and rationale.

### 7.1 Tool inventory

New tools sourced from BioContainers (pinned-hash):

| Tool | Bioconda recipe | Container URI |
|---|---|---|
| freebayes 1.3.10 | `freebayes=1.3.10` | `quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0` *(landed M1; PR #70 corrected the registry prefix from `biocontainers/` to `quay.io/biocontainers/`)* |
| LoFreq 2.1.5 | `lofreq=2.1.5` | `quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4` *(landed M2)* |
| VarDict-Java 1.8.3 | `vardict-java=1.8.3` (bundles `teststrandbias.R` + `var2vcf_valid.pl`) | `quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0` *(landed M2 via COMB-418)* |
| bedtools 2.31.1 | `bedtools=2.31.1` | `quay.io/biocontainers/bedtools:2.31.1--h13024bc_3` *(landed M2 via COMB-418 ; powers `BEDTOOLS_MAKEWINDOWS`)* |
| htslib 1.21 (bgzip) | `htslib=1.21` | reused via `community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11` *(landed M2 via COMB-418 ; powers `BGZIP` for callers whose own image lacks bgzip, currently `vardict-java`)* |
| GATK 4.6.1.0 | `gatk4=4.6.1.0` | `quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0` *(landed M2 via COMB-419 / PR #70 ; powers both `GATK4_HAPLOTYPECALLER` and `GATK4_CREATESEQUENCEDICTIONARY`)* |

Consumed directly (upstream-published, not BioContainers):
- **DeepVariant** ; `google/deepvariant:1.6.1` (official upstream image). Large (~7 GB) but production-grade and self-contained. Do not rebuild ; Google's image is the canonical distribution and bundles trained models. `TODO: check the licenses for the trained models and if we can use them`

Reused (no action):
- `bcftools:1.21` ; for `BCFTOOLS` (caller), `BCFTOOLS_NORM`, `BCFTOOLS_STATS`, and the `samtools faidx` module (samtools is bundled in this image).
- `samtools` / `tabix` ; bundled in existing images.
- polars + pysam ; reuse `pip_pyyaml_upsetplot_numpy_pruned` for `VCF_TO_PARQUET`, `COUNT_VARIANT`, metrics YAML.

**Rationale for `vardict-java` over perl VarDict:** Higher throughput.

**Container risks:** A BioContainers tag may be retracted upstream if the recipe is unmaintained (low probability — pinned hashes survive). Mitigation: if it ever happens, fall back to a Wave rebuild of the same Bioconda recipe (§7.2 step 3).

### 7.2 Creating and using containers

This subsection is the general guide for "I need a container for tool X". The decision steps below are ordered from easiest to most difficult; only fall through if the previous option doesn't fit. **For variant-calling tools (v1), step 2 is the canonical choice — BioContainers pinned-hash images.** Wave (step 3) is the fallback if a BioContainers image is unavailable or insufficient.

**1. Reuse an image that already exists in the repo.** Grep for the tool first ; if `bcftools`, `samtools`, `tabix`, `pysam`, `cyvcf2`, `polars`, `pyyaml`, or `multiqc` is what we need, an existing module is already pulling an image that has it. For variant calling specifically:
   - Anything that calls `bcftools` / `samtools` / `tabix` ; reuse `community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11` ([modules/dna_processing.nf:82](modules/dna_processing.nf#L82)).
   - Custom Python (pysam, cyvcf2, polars, pyyaml) ; reuse `community.wave.seqera.io/library/pip_pyyaml_upsetplot_numpy_pruned:3f7a3b0ec65284cc` ([modules/shared.nf:3](modules/shared.nf#L3)). This is the image referenced in [README.md:131](README.md#L131); its Wave build record is at [wave.seqera.io/view/builds/bd-3f7a3b0ec65284cc_1](https://wave.seqera.io/view/builds/bd-3f7a3b0ec65284cc_1).

**2. Use a BioContainers pinned-hash image.** This is the v1 default for variant-calling tools (freebayes, lofreq, vardict-java, gatk4):
   - Search [BioContainers on quay.io](https://quay.io/organization/biocontainers) (or the [BioContainers registry](https://biocontainers.pro/)) for the tool and version.
   - Copy the pinned-hash tag of the form `<version>--<build>` (e.g. `1.3.10--hbefcdb2_0` for freebayes 1.3.10).
   - Paste `quay.io/biocontainers/<tool>:<version>--<build>` into the module's `container "..."` line. Example: `quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0` ([modules/freebayes/main.nf:4](modules/freebayes/main.nf#L4)). **Always include the `quay.io/` registry prefix** ; the bare `biocontainers/<tool>:<tag>` form resolves against Docker Hub and fails to pull (caught in COMB-419 / PR #70; freebayes was originally pasted in without the prefix and `docker pull` returned `not found`).
   - Never use a floating tag (e.g. `quay.io/biocontainers/freebayes:latest` or `:1.3.10` without the build suffix). The `--<build>` hash is the content-addressed pin.

**3. Wave-build a single bioconda package.** Fallback for tools without a usable BioContainers image, or when a custom layering is needed:
   - Open the [Wave container builder](https://seqera.io/containers/).
   - Select Conda packages, paste e.g. `bioconda::freebayes=1.3.10`.
   - Build, copy the resulting URI of the form `community.wave.seqera.io/library/freebayes:1.3.10--<hash>`, paste into the module's `container "..."` line.
   - The hash pins the build content-addressably ; IMPORTANT for reproducibility: never use a floating tag!!!

**4. Wave-build a multi-tool image.** Same flow as (3) but with multiple packages in one request. The published name will be hash-suffixed (`community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c` is an existing example in [modules/dna_processing.nf:131](modules/dna_processing.nf#L131)). Use this only when the tools are coupled inside a single process body (`bwa | samtools sort` style); otherwise prefer separate processes with separate images so each step is independently cacheable and resumable.

**5. Use an upstream-published image directly.** When the upstream project ships a canonical image with bundled assets, reference it directly instead of rebuilding. This is the case for **DeepVariant** ; `google/deepvariant:1.6.1` ships pre-trained models that a BioContainers / Wave rebuild would not. Pin a specific version (`:1.6.1`), never `:latest`. Document the choice in the module body.

**6. Build a custom image (Python helpers, R, shell tooling).** Same Wave UI flow as (3) but submit either a Conda env spec, a `pip`-packages list, or a `Dockerfile`. The output is again a `community.wave.seqera.io/library/...:<hash>` URI to paste into the module.

**Reproducibility gap (proposed action item).** The recipes that produced the *custom* Wave images currently live only as Wave build records, not in this repo. If a Wave build record is ever lost or the source recipe drifts, the custom image is unreconstructible from the repo alone. As part of Milestone 1, add a `containers/` directory containing the source for each custom image (`environment.yml` / `requirements.txt` and/or `Dockerfile`) so the image can be rebuilt deterministically. The pinned `community.wave.seqera.io/...` and `quay.io/biocontainers/...` URIs in module bodies will remain the source of truth.

**Update procedure for an existing image (e.g. `freebayes 1.3.10` → `1.3.11`, or `bcftools 1.21` → `1.22`).**
1. For BioContainers tools: locate the new pinned-hash tag on quay.io and substitute. For Wave tools: rebuild via the UI or `wave` CLI.
2. Open one PR per tool updating the URI(s); avoid mass-updating multiple tools in one PR so blast radius stays small if the new image regresses.
3. Run the §8.2 sanity test before merge.
4. Pin the new content hash; never replace a hash with a floating tag during the bump.

**Why BioContainers for v1 variant calling (and not all-Wave)?** The dual-URI `quay.io/biocontainers/...` + `https://depot.galaxyproject.org/singularity/...` pattern that nf-core uses is upstream-maintained, content-addressed, and avoids the Wave-build round-trip on every tool addition. For the v1 variant-calling work (freebayes, plus the four M2/M3 callers) BioContainers ships everything we need at the right pins, so we adopt it here. The pre-existing pipeline modules stay on Wave (no churn); revisiting whether to migrate those to BioContainers is a future refactor.

### 7.3 Tool maintainability risk

As we discussed in the risk registry, each third-party tool carries a maintainability risk. The table below records comments, the consequence of it going stale, and a potential mitigation plan. `TODO: re-review this table at every major uniflow release or maybe more often?`

**Summary of the risk:**
- The highest single-tool risk is **LoFreq**. The architecture keeps each caller structurally decoupled via `params.variant_callers`, so LoFreq can be retired without touching any other module.
- Two of the six callers are semi-active or in maintenance mode (`freebayes`, `VarDict-Java`). Our fall back scenarios can cover for these two.
- The two newest additions (`DeepVariant`, `GATK_HC`) are both low-risk and backed by Google, Broad Institute.
- Normalization/QC infrastructure (`bcftools`) is low-risk.

| Tool | Version pin | Last upstream release (approx.) | Upstream status | Risk | Exit plan if abandoned |
|---|---|---|---|---|---|
| `bcftools` / `samtools` / `htslib` | 1.21 | 2024 | **Active** ; Sanger Institute; regular releases, large contributor base | **Low** | No realistic replacement needed; core community dependency. |
| `freebayes` | 1.3.10 | 2024 | **Semi-active** ; slow, seems there are years between tags | **Moderate/High** | Lean more heavily on DeepVariant or GATK HaplotypeCaller; both are already in the default caller set. |
| `LoFreq` (LoFreq2) | 2.1.5 | ~2018 | **Stale** ; no releases since 2018; LoFreq3 exists in the author's repo but has no public release | **Very High** | Primary replacement path: lower VarDict's `-f` threshold to ~0.01 to cover the low-frequency niche. Longer term: DeepVariant can replace LoFreq's "low-VAF sensitivity" role in many contexts, or adopt LoFreq3 in a world where we might get a new release!. If the bioconda recipe stops building, drop LoFreq from the default `params.variant_callers` rather than blocking a release. Or completely drop LoFreq from the bundle. |
| `VarDict-Java` | 1.8.3 | ~2020 | **Stale** ; No new releases will be made (same for `VarDict Perl`) | **Very High** | If the maintainability doesn't change, fall back to GATK Mutect2 germline mode (already containerized via the `gatkhc` BioContainers image). |
| `DeepVariant` | 1.6.1 | 2024 | **Active** ; Google Health; frequent releases, production-grade models for Illumina/PacBio/ONT | **Low** | Drop from the default caller set if for example license changes and becomes more restrictive. |
| `GATK HaplotypeCaller` | 4.6.1.0 | 2024 | **Active** ; Broad Institute; industry-standard germline caller, regular quarterly releases | **Low** | Drop from the default caller set, or pin an older GATK version. HaplotypeCaller has been stable for years. |
| `pysam` / `cyvcf2` / `polars` / `multiqc` | current | 2024 | **Active** | **Low** | Standard Python/bioinformatics dependencies; replaceable with other Rust based alternatives. |

---

## 8. Testing Strategy

Point of concern and question: How reliable is creating dummy fixture for Nextflow similar to Python. If it is possible, we should have this as part of the CI or pipelines on the BitBucket.
`TODO: investigate how we can create dummy fixtures for Nextflow as part of the CI`

### 8.1 Unit / script-level
- `bin/vcf_to_parquet.py` ; A dummy minimal VCF (SNV, multi-allelic, missing AD, missing AF) to parquet rows.
- `bin/variants_to_metrics_parquet.py` ; A dummy parquet to expected YAML.
`TODO: expand testing strategy`

### 8.2 Integration (Nextflow)
- **Stub-run topology check:** `nextflow run main.nf … -stub-run` exercises the entire DAG end-to-end without running real tool bodies. Every process declares a `stub:` block that `touch`es files matching its declared output globs, so channel topology, joins, optional-output binding, and the per-sample `VARIANT_CALLING` filter are validated without container pulls or compute. This is the gate that catches the class of bugs `-preview` does not (empty-output silent failures, wrong stub filenames, schema/meta-shape mismatches like the `nf-schema` empty-boolean materialisation).
- **Sanity check:** `nextflow run main.nf -profile docker,test --variant_callers freebayes …` asserts:
  - `samples/*/DNA/variants/freebayes/*.vcf.gz` exists.
  - `bcftools view -h` parses cleanly.
  - `bcftools view -v snps $vcf | bcftools view -H | wc -l` equals `bcftools view -H $vcf | wc -l` (SNV-only invariant).
  - `PILEUP_READS` is absent from the resumed `.nextflow.log` (legacy-removal invariant).
- **Full six-caller (if we decide to use them all):** all six VCFs + six parquets per sample; MultiQC shows six per-caller `bcftools stats` sections.
- **Per-cell (if we decide to produce per cell VCF):** per-caller `variants_count.h5ad` non-empty on fixture.

### 8.3 CI
[bitbucket-pipelines.yml](bitbucket-pipelines.yml) runs the §8.2 stub-run on every PR (replaces the prior `-preview` step). It exercises the full DAG including the per-sample `run_variant_calling` filter against a fixture samplesheet covering `true`/`false`/empty cells and an RNA row. The real-tool sanity check (containers, fixture FASTQ) is reserved for a later, heavier CI step. DeepVariant and GATK are **not** in the CI default caller set ; their images are too heavy. The eventual real-tool sanity will run `--variant_callers freebayes,bcftools`; full six-caller coverage (if we decide to use all 6) can run on another test.

---

## 9. Implementation Plan and milestones

### Milestone 1 ; Skeleton + freebayes + legacy code removal
Deliverables:
- [x] `params.variant_callers` (schema).
- [x] `params.run_variant_calling` (schema; default `true`) and samplesheet `run_variant_calling` column ([assets/schema_input.json](assets/schema_input.json)).
- [x] Filter at the `DNA_PROCESSING to VARIANT_CALLING` level by applying the logic (per-sample override, falling back to the global param). Note: `nf-schema` materialises an empty optional-boolean cell as `[]` (not `null`/`false`); the filter therefore checks `instanceof Boolean` rather than `!= null` to keep "empty cell → fallback to `params`" working.
- [x] Per-tool modules: [modules/samtools/faidx/main.nf](modules/samtools/faidx/main.nf) (`SAMTOOLS_FAIDX`, COMB-402), [modules/bcftools/norm/main.nf](modules/bcftools/norm/main.nf) (`BCFTOOLS_NORM`, COMB-403), [modules/bcftools/stats/main.nf](modules/bcftools/stats/main.nf) (`BCFTOOLS_STATS`, COMB-405), [modules/tabix/main.nf](modules/tabix/main.nf) (`TABIX`), [modules/freebayes/main.nf](modules/freebayes/main.nf) (`FREEBAYES`, COMB-404). *(`FREEBAYES` now points at the BioContainers image `quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0`; QUAL/DP filtering inside `BCFTOOLS_NORM` is deferred to the schema-params task COMB-410.)*
- [x] [subworkflows/variant_calling.nf](subworkflows/variant_calling.nf) single-caller path with `(meta, caller)` channel shape.
- [x] Wiring into `DNA_PROCESSING`: invoke `VARIANT_CALLING` (COMB-406, on branch `feature/COMB-406-connect-variant_calling-workflo3`, not yet merged).
- [x] **Disable `COUNT_VARIANT` and its emissions (`snp_h5ad`, `snp_dev`, `snp_metrics`) in `DNA_PROCESSING`** (COMB-407). Restored in M4 once the merged-VCF parquet catalog (§4.7) exists; in M1 the per-cell attribution layer is intentionally absent because pileup removal orphans its sole input. Process body and `bin/count_variant.py` retained for M4 rewiring.
- [x] **Remove `PILEUP_READS` from `modules/dna_processing.nf` and all references** (COMB-407). Legacy `snps.csv` publish path removed from [main.nf](main.nf), [subworkflows/dna_processing.nf](subworkflows/dna_processing.nf), [subworkflows/cross_processing.nf](subworkflows/cross_processing.nf), and the [README](README.md) example tree.
- [x] Update [main.nf](main.nf) with publish paths for `samples/{sample}/DNA/variants/{caller}/` (COMB-408). `DNA_PROCESSING` emits `variant_outputs` (mixed `vcf` + `tbi` + `stats`, grouped by `(meta, caller)`); main.nf publishes via `path { meta, caller, files -> "samples/${meta.sample_name}/DNA/variants/${caller}/" }`.
- [x] Schema validation: `variant_callers: []` rejected with a message pointing to `run_variant_calling` ([nextflow_schema.json:99-105](nextflow_schema.json#L99-L105)). Remaining schema work for `variant_min_qual`, `variant_min_dp`, `raw_vcf` is COMB-410.
- [x] FreeBayes container: `quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0` (COMB-409; registry prefix corrected to `quay.io/biocontainers/` in COMB-419 / PR #70). Version pinned at 1.3.10 (vs. originally planned 1.3.7); switched from a Wave-built URI to a BioContainers pinned-hash image — see §7.
- [x] Sanity test runnable on fixture (`--variant_callers freebayes`) covering the four assertions in §8.2.

Acceptance criteria:
- sanity test passes. Repo has no references to the removed legacy outputs. MultiQC shows one bcftools sections.
- A DNA sample with `run_variant_calling: false` (samplesheet) produces no `samples/{sample}/DNA/variants/...` output and produces no `CALL_*` / `BCFTOOLS_NORM` / `COUNT_VARIANT` tasks in `.nextflow.log`. Other DNA outputs (alignment, amplicon counts, QC) still produced.
- `params.run_variant_calling = false` skips variant calling entirely for every DNA sample regardless of the samplesheet column.
- `variant_callers: []` is rejected at schema validation with a message pointing to `run_variant_calling`.

Validation status (2026-05-06, COMB-407 stub-run on a 4-row fixture: S1_DNA `true`, S2_DNA `false`, S3_DNA empty, S4_RNA empty):
- [x] Per-sample filter behaviour: `FREEBAYES 2 of 2` (S1_DNA + S3_DNA via `params` fallback); S2_DNA correctly excluded. `nf-schema` empty-cell `instanceof Boolean` check holds.
- [x] Legacy-removal invariant: no `PILEUP_READS` / `COUNT_VARIANT` tasks in the trace; no references to `_snps.csv` / `filtered_vcf_csv` / `snp_h5ad` / `snp_dev` / `snp_metrics` outside `docs/`.
- [x] `VARIANT_CALLING` topology: `SAMTOOLS_FAIDX` → `FREEBAYES` → `BCFTOOLS_NORM` → `TABIX` + `BCFTOOLS_STATS` wires through end-to-end.
- [x] Publish-layer per-sample isolation (COMB-408 stub-run): S1_DNA and S3_DNA have `samples/{sample}/DNA/variants/freebayes/{vcf.gz, tbi, bcftools_stats.txt}`; S2_DNA's `DNA/` dir contains only alignment outputs, no `variants/` subdirectory.
- [x] `params.run_variant_calling = false` global-skip path: not yet exercised; needs a second stub-run with `--run_variant_calling false`.

Real-data E2E (2026-05-08, COMB-412 merged in PR #65):
- [x] Full pipeline run on actual data (single-sample fixture beyond stub-run); `--variant_callers freebayes` produces non-empty VCFs and bcftools stats without runtime errors. Confirms the M1 skeleton plus FreeBayes container (COMB-409) survive a real BAM and reference FASTA, not just a `touch`-stubbed DAG.


### Milestone 2 ; Remaining CPU callers + parquet

Scope: add the four remaining CPU-only callers (LoFreq, bcftools, VarDict, GATK_HC) alongside the M1 freebayes path, land the BED fallback that VarDict + GATK_HC depend on, ship the parquet emitter, and wire per-caller bcftools stats into MultiQC. M2 does **not** add the per-cell attribution layer (M4) or DeepVariant (M3).

Deliverables (per-tool modules, each at `modules/<tool>/main.nf` ; one ticket per module ; "standalone per-tool module" means each caller is its own process file, not a multi-tool umbrella). Sequencing matters: items below are listed in execution order. The MultiQC + metrics-path move comes first because once it lands, every subsequent caller's `bcftools stats` output is auto-picked-up with no further wiring.

- [x] **MultiQC wiring + bcftools stats publish-path move (FIRST)** ; per §4.5 "MultiQC wiring + publish-path move" and §5. Splits `BCFTOOLS_STATS` output out of `variant_outputs` (which becomes vcf+tbi only), adds a dedicated `variant_stats` emit on `DNA_PROCESSING`, mixes that channel into `ch_multiqc_files` in [main.nf](main.nf), and publishes stats files at the new flat location `samples/{sample}/metrics/variants/`. Touches [main.nf](main.nf) and [subworkflows/dna_processing.nf](subworkflows/dna_processing.nf). Refactor of the M1 freebayes path: the freebayes stats file moves out of `samples/{sample}/DNA/variants/freebayes/` (M1 alpha; no external consumers, no compat shim). Once this lands, every M2 caller's stats file is mixed in automatically because `BCFTOOLS_STATS` already emits `(meta, caller, stats)` shaped channels for any caller normalized through `BCFTOOLS_NORM`. Validated via stub-run on the 4-row fixture (2026-05-08): S1_DNA + S3_DNA each publish `samples/{sample}/metrics/variants/{sample}.freebayes.bcftools_stats.txt`; S2_DNA (`run_variant_calling=false`) has no `metrics/variants/` dir; `BCFTOOLS_STATS 2 of 2` confirms the channel reaches MultiQC for the two enabled samples only.
- [x] **`LOFREQ_CALL`** ([modules/lofreq/call/main.nf](modules/lofreq/call/main.nf)) ; `lofreq call-parallel --pp-threads $task.cpus -f $fasta -o ${prefix}.vcf.gz $bam` per §4.2. No `--call-indels`, no `lofreq indelqual` prereq. `cpu_low_mem` label. Container: `quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4` (§7.1). Output tuple: `(meta, "lofreq", path("${meta.sample_name}.lofreq.vcf.gz"))`. Per the "LoFreq notes" in §4.2, the lofreq VCF is missing `##contig=` header lines, so the stream is reheadered via `GATK4_UPDATEVCFSEQUENCEDICTIONARY` ([modules/gatk4/updatevcfsequencedictionary/main.nf](modules/gatk4/updatevcfsequencedictionary/main.nf), `gatk IndexFeatureFile` + `UpdateVCFSequenceDictionary --source-dictionary <bam> --replace true`) before entering `caller_vcfs`. The reheader module reuses the GATK 4.6.1.0 BioContainers image already pinned for `CALL_GATKHC`; no new container. The lofreq stream is wired in [subworkflows/variant_calling.nf](subworkflows/variant_calling.nf) by joining `LOFREQ_CALL.out.vcf` with `bam_with_index` on meta. As a side effect, `TABIX` was switched to `tabix --force` so it tolerates a pre-existing `.tbi` in the work dir. Module path moved from `modules/lofreq/main.nf` to `modules/lofreq/call/main.nf` to match the per-subcommand layout used by other multi-step tool families (`bcftools/`, `gatk4/`).
- [ ] **`CALL_BCFTOOLS`** ; `bcftools mpileup -a AD,DP,SP --skip-indels -f $fasta $bam | bcftools call -mv --ploidy 2 -Oz -o out.vcf.gz` per §4.2. `--skip-indels` at the `mpileup` level makes the `bcftools view -v snps` step in `BCFTOOLS_NORM` redundant for this caller, which is fine ; the duplication keeps the module pure. Reuses the existing `bcftools:1.21` Wave image (§7); no new container. `small_job` or `cpu_low_mem` label. Output tuple: `(meta, "bcftools", path("${meta.sample_name}.bcftools.vcf.gz"))`.
- [x] **`CALL_VARDICT`** (COMB-418) ; `vardict-java -f 0.2 -r 2 -c 1 -S 2 -E 3 -g 4 -th $task.cpus -G $fasta -b $bam $amplicon_bed | teststrandbias.R | var2vcf_valid.pl -A -f 0.2 -N $meta.sample_name > ${prefix}.vcf` per §4.2 ([modules/vardict/main.nf](modules/vardict/main.nf)). `vardict-java`, not perl VarDict. `var2vcf_valid.pl -N` takes the sample name from `meta.sample_name`. `cpu_low_mem` label. Container: `quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0`. Consumes `amplicon_bed_ch` (§4.3.3). The image does not bundle `bgzip`, so the module emits an uncompressed `.vcf` and the subworkflow runs `BGZIP` ([modules/tabix/bgzip/main.nf](modules/tabix/bgzip/main.nf)) before joining `caller_vcfs`. Output tuple after BGZIP: `(meta, "vardict", path("${meta.sample_name}.vardict.vcf.gz"))`.
- [x] **`CALL_GATKHC`** (COMB-419, PR #70) ; `gatk --java-options "-Xmx8g" HaplotypeCaller --reference $fasta --input $bam --output out.vcf.gz --tmp-dir . --output-mode EMIT_VARIANTS_ONLY --sample-ploidy 2` ([modules/gatk4/haplotypecaller/main.nf](modules/gatk4/haplotypecaller/main.nf)). **`-Xmx8g` hardcoded inside the module body**, not threaded from `task.memory`; rationale in §4.2 GATK_HC notes. `cpu_low_mem` label. Container: `quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0`. Consumes `GATK4_CREATESEQUENCEDICTIONARY.out.dict`. Output tuple: `(meta, "gatkhc", path("${meta.sample_name}.gatkhc.vcf.gz"))`. **`-L $amplicon_bed` is still not wired** ; `BEDTOOLS_MAKEWINDOWS` landed in COMB-418 but the GATK_HC wiring was kept out of scope. Add the BED flag (and `amplicon_bed_ch` as an input on `CALL_GATKHC`) in a follow-up; extend `bed_consumers` in [subworkflows/variant_calling.nf](subworkflows/variant_calling.nf) to include `'gatkhc'` so the BED channel is materialised whenever GATK_HC is requested.
- [x] **`GATK4_CREATESEQUENCEDICTIONARY`** (COMB-419, PR #70) ; `gatk --java-options "-Xmx6g" CreateSequenceDictionary --REFERENCE ${fasta} --URI ${fasta} --TMP_DIR .` producing `${fasta.baseName}.dict` ([modules/gatk4/createsequencedictionary/main.nf](modules/gatk4/createsequencedictionary/main.nf)). `small_job` label. Reuses the GATK 4.6.1.0 BioContainers image already pinned for `CALL_GATKHC`; no new container. Runs once inside `VARIANT_CALLING`, gated on `gatkhc in callers`.
- [x] **`BEDTOOLS_MAKEWINDOWS`** (COMB-418) ; full spec at §4.3.3. Conditional on `params.amplicon_bed == null`; consumed by `CALL_VARDICT` (and later `CALL_GATKHC` and `CALL_DEEPVARIANT`). Lives at [modules/bedtools/makewindows/main.nf](modules/bedtools/makewindows/main.nf). Implementation switched from the original `awk` sketch to `bedtools makewindows -g <fai> -n 1`; **adds one new container** (`quay.io/biocontainers/bedtools:2.31.1--h13024bc_3`), deviating from the original "no new container" plan ; rationale in §4.3.3.
- [ ] **`VCF_TO_PARQUET`** + **`bin/vcf_to_parquet.py`** ; full schema at §4.7. Reuses `community.wave.seqera.io/library/pip_pyyaml_upsetplot_numpy_pruned:3f7a3b0ec65284cc` (verify `cyvcf2`/`pysam` is present in the image; if absent, layer in via Wave rebuild ; one of those two should be ; see §7.2). Process consumes `(meta, caller, vcf)` from `BCFTOOLS_NORM`; emits `(meta, caller, path("${meta.sample_name}.${caller}.variants.parquet"))`. Empty VCF → empty-but-schema-correct parquet (zero rows, 12 columns); not an error.
- [ ] **Schema params** ; `params.variant_min_qual` (object<string,int>), `params.variant_min_dp` (integer), `params.raw_vcf` (boolean) per §6.1. Plumbed into `BCFTOOLS_NORM` (`variant_min_qual` per-caller, `variant_min_dp` site-level; both now active in the §4.4 command body, not deferred), and into the publish layer (`raw_vcf=true` publishes pre-norm VCFs at `samples/{sample}/DNA/variants/{caller}/{sample}.{caller}.raw.vcf.gz` per §5).
- [x] **Containers** ; VarDict-Java 1.8.3 pinned (`quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0`, COMB-418); bedtools 2.31.1 added for `BEDTOOLS_MAKEWINDOWS` (`quay.io/biocontainers/bedtools:2.31.1--h13024bc_3`, COMB-418); LoFreq 2.1.5 and GATK 4.6.1.0 already pinned (§7.1). The §7 default is BioContainers pinned-hash images; if a usable hash is unavailable or fails on the fixture, fall through to a Wave build for that single tool (§7.2 step 3). The bcftools-only `CALL_BCFTOOLS` reuses the existing image and does not need a new container.
- [ ] **Real-tool sanity** ; the four §8.2 assertions, plus the MultiQC bcftools section count. Carried over from M1's deferred item; running the real-tool fixture against the M2 caller set (`--variant_callers freebayes,lofreq,bcftools,vardict,gatkhc`) covers both. M1 has been validated end-to-end on real data via COMB-412, so this is now strictly an additive check for the four new callers.

Open decisions inside M2 (resolve at PR time, not now):
- LoFreq parallel sharding granularity on amplicon FASTAs (`call-parallel` shards by contig by default; with many short amplicon contigs this should saturate `cpu_low_mem` ; revisit if profiling shows otherwise).
- Whether `CALL_BCFTOOLS` lives at `cpu_low_mem` or `small_job`. Default to `cpu_low_mem` to match the other callers; downsize if the mpileup `|` call pipe is consistently underutilising threads on the fixture.

Acceptance criteria:
- All five CPU-caller VCFs (`freebayes`, `lofreq`, `bcftools`, `vardict`, `gatkhc`) and five matching parquets produced per sample on the real-data fixture. VarDict and GATK_HC each produce non-empty output (the indel-heavy edge cases get filtered in `BCFTOOLS_NORM`, but at least one SNV record must remain).
- MultiQC report shows five `bcftools stats` sections per sample, one per caller, auto-discovered via the `<sample>.<caller>.bcftools_stats.txt` filename pattern.
- `params.amplicon_bed = null` runs `BEDTOOLS_MAKEWINDOWS` once and feeds `CALL_VARDICT` from the derived BED (and `CALL_GATKHC` once that wiring lands ; see §10 GATK_HC entry). `params.amplicon_bed = <file>` skips `BEDTOOLS_MAKEWINDOWS` (verifiable via `.nextflow.log` task absence). When no BED-consuming caller is in `params.variant_callers`, `BEDTOOLS_MAKEWINDOWS` is not invoked even if `params.amplicon_bed` is null.
- `raw_vcf = true` publishes `{sample}.{caller}.raw.vcf.gz` alongside the normalized VCF for each of the five callers; `raw_vcf = false` omits them. The normalized VCF is published either way.
- `variant_min_qual` / `variant_min_dp` are observable in `BCFTOOLS_NORM`'s VCF: a record with `QUAL < threshold` or `INFO/DP < params.variant_min_dp` does not appear in the normalized output. (Spot-check on the fixture; not a full filter unit test.)
- `--variant_callers freebayes` (M1's path) still produces exactly the M1 outputs ; M2 is additive, not regressive.
- Per-sample `run_variant_calling: false` and global `params.run_variant_calling = false` continue to skip the whole subworkflow including all M2 additions.

Validation status (2026-05-11, COMB-419 / PR #70 — GATK HaplotypeCaller + sequence dictionary landed):
- [x] `GATK4_CREATESEQUENCEDICTIONARY` produces `${fasta.baseName}.dict` on the amplicon FASTA fixture; runs once per pipeline invocation when `gatkhc` is in `params.variant_callers`.
- [x] `GATK4_HAPLOTYPECALLER` produces `${meta.sample_name}.gatkhc.vcf.gz` on the real-data fixture and feeds `BCFTOOLS_NORM` → `TABIX` + `BCFTOOLS_STATS` cleanly. **One-time fix required:** GATK rejected the BAM with "the sample list cannot be null or empty" because `MINIMAP2`'s output had no `@RG`; resolved by adding a read group at alignment time (see commit history on PR #70). All four §8.2 sanity assertions hold for the `gatkhc` caller stream.
- [x] Schema default updated to `params.variant_callers = "freebayes,bcftools,gatkhc"` ([nextflow_schema.json](nextflow_schema.json)); the GATK HC stream is now exercised on every default-config run, not just opt-in invocations.
- [ ] BED restriction (`-L $amplicon_bed`) on `CALL_GATKHC` still not wired. `BEDTOOLS_MAKEWINDOWS` itself landed in COMB-418 ([modules/bedtools/makewindows/main.nf](modules/bedtools/makewindows/main.nf), `bedtools makewindows -g <fai> -n 1`) and is wired into `CALL_VARDICT`, but `CALL_GATKHC` was kept out of scope. Follow-up: add `path amplicon_bed` to `GATK4_HAPLOTYPECALLER`, pass `-L ${amplicon_bed}`, and extend `bed_consumers` in [subworkflows/variant_calling.nf](subworkflows/variant_calling.nf) to include `'gatkhc'`.

### Milestone 3 ; DeepVariant
Deliverables:
- `CALL_DEEPVARIANT` consuming `google/deepvariant:1.6.1`.
- `params.deepvariant_model_type` plumbed (default `WES`).
- Resource label review (DeepVariant is heavier than CPU callers).

Acceptance criteria:
- Six VCFs + six parquets per sample on the fixture. MultiQC shows six bcftools sections.

### Milestone 4 ; Cross-caller merge + per-cell attribution (SNV)

Deliverables:

**Cross-caller VCF surface (in progress).**
- [ ] **Merge per-caller VCFs into a single sample-level VCF (COMB-439).** Harmonize → reheader → merge → annotate path landed in [subworkflows/variant_calling.nf](subworkflows/variant_calling.nf) ; design in [docs/adr/variant_merging/variant_merging.md](docs/adr/variant_merging/variant_merging.md). Lofreq is dormant in the default `variant_callers` (its site-only output lacks FORMAT/AD,DP — no clean bcftools-native lift exists, deferred), but the template + merge plumbing is retained for revival. Per-caller `BCFTOOLS_FILL_TAGS` ([modules/bcftools/fill_tags/main.nf](modules/bcftools/fill_tags/main.nf)) only copies `QUAL` into `INFO/VARQUAL` so per-caller site quality survives the merge. The consensus VCF retains per-caller `FORMAT/{VAF,AD,DP,GT}` columns and carries cross-caller summaries in INFO (`VAF_{min,max}`, `NCALLERS`, `DP_{min,median,max}`, `AD_REF_{min,median,max}`, `AD_ALT_{min,median,max}`) — all synthesised post-merge via a `bcftools +fill-tags` pipe in `BCFTOOLS_ANNOTATE_CONSENSUS` (no custom scripts). New modules: [modules/bcftools/fill_tags/main.nf](modules/bcftools/fill_tags/main.nf), [modules/bcftools/reheader_sample/main.nf](modules/bcftools/reheader_sample/main.nf), [modules/bcftools/merge/main.nf](modules/bcftools/merge/main.nf), [modules/bcftools/annotate_consensus/main.nf](modules/bcftools/annotate_consensus/main.nf), [modules/template_vcf/main.nf](modules/template_vcf/main.nf) (header-only VCF used to attach a sample column to lofreq's site-only output via `bcftools merge`). [modules/bcftools/norm/main.nf](modules/bcftools/norm/main.nf) switched to `task.ext.args` (adds `--rm-dup all`). New [conf/modules.config](conf/modules.config) introduced (sarek-style `withName` + `ext.args`); only new/touched modules adopt the pattern. Published to `samples/{sample}/DNA/variants/merged/`.
- [ ] **Summarize and compare per-caller VCF outputs (COMB-438).** Shipped together with COMB-439 as `CALLER_COMPARE_QC` ([modules/qc/caller_compare/main.nf](modules/qc/caller_compare/main.nf)), a Nextflow wrap of [bin/compare_caller_runs.sh](bin/compare_caller_runs.sh) (copied verbatim from [docs/adr/lofreq/compare_caller_runs.sh](docs/adr/lofreq/compare_caller_runs.sh); flagged for decomposition into smaller pieces as a follow-up, the single-script form is too large to maintain). Publishes `pairwise_jaccard.tsv`, `af_pairwise_stats.tsv`, `summary.txt`, etc. under `samples/{sample}/qc/variants/`.

**Per-cell attribution + filtering surface.**
- [x] **Schema params `variant_min_dp` / `variant_min_vaf` (COMB-433 + COMB-436).** Per §6.1. `variant_min_dp` plumbed into the new `BCFTOOLS_FILTER` (§4.4.1) as a site-level `INFO/DP` floor downstream of `BCFTOOLS_NORM`; per-cell DP gating in `count_variant.py` removed entirely (the parquet's per-caller `FORMAT/DP` is a caller-side report, not a per-cell count). `variant_min_qual` dropped from scope — per-caller QUAL scales differ too much for a single defensible cutoff without truth data, and per-caller cutoffs are deferred until a benchmarking fixture exists. `variant_min_vaf` added in its place as a per-cell hard floor inside `filter_snps()`, intended to suppress low-VAF sequencing/PCR/mapping noise without touching the homo/hete classification.
- [x] **`VCF_TO_COUNTS_CSV` (COMB-432) — SUPERSEDED.** Collapsed into `VCF_TO_PARQUET` (§4.6, §4.7). The long-format parquet is a strict superset of the originally specified CSV schema and per-caller QUAL comes directly from `INFO/VARQUAL[idx]`, so the CSV intermediate adds nothing.
- [x] **`COUNT_VARIANT` per-caller fan-out (COMB-434).** Per §4.8. Reuses the existing `COUNT_VARIANT` body and [bin/count_variant.py](bin/count_variant.py); rewired to consume the merged-VCF-derived parquet (filtered by `caller`) instead of the removed pileup CSV. Runs once per `(meta, caller)`.
- [x] **Refactor hardcoded values in `count_variant.py` (COMB-436).** `MIN_BASE_QUALITY=30` (a misnomer — it was filtering the parquet's VCF QUAL, not base quality) and `MIN_DP=200` both removed; the QUAL filter is gone (per `variant_min_qual` decision above) and the DP filter migrated upstream to `BCFTOOLS_FILTER`. The actual base-quality knob inside `pileup_truncated` (`min_base_quality=13`) also removed — set to `0` because the pileup is an attribution step and per-base quality belongs to the callers (§4.2 "Base-quality filtering"). `MIN_TOTAL_READS=4` and `FRACTION_READS_HOMO=0.8` retained as module constants with a docstring note marking them as candidate future params (per-cell genotype-calling knobs, not variant-quality gates — promote once a truth set justifies tuning).
- [x] **Publish per-caller `variants_count.h5ad` (COMB-435).** Published under `samples/{sample}/DNA/variants/{caller}/per_cell/` per §5 via the new `variant_per_cell` channel emitted from `DNA` and routed in [main.nf](main.nf).
- [x] **Investigate additional counts-CSV columns (COMB-437) — SUPERSEDED.** Resolved by the parquet schema in §4.7: per-caller `FORMAT/AD` (`ad_ref`/`ad_alt`), per-caller VAF (`af`), and `FILTER` are all first-class columns. `caller_specific` is a placeholder string in v1, to revisit if/when per-caller INFO is preserved through `BCFTOOLS_MERGE`.

**Parquet.**
<<<<<<< HEAD
- [x] **`VCF_TO_PARQUET` + `bin/vcf_to_parquet.py` (COMB-422).** Per §4.7. Long-format parquet sourced from the merged consensus VCF, one row per (sample, caller, variant), 12 columns. Per-caller QUAL via `INFO/VARQUAL[idx]`. Doubles as the catalog source for `COUNT_VARIANT`.
=======
- [x] **`VCF_TO_PARQUET` + `bin/vcf_to_parquet.py` (COMB-422).** Per §4.7. Long-format parquet sourced from the merged consensus VCF, one row per (sample, caller, variant), 13 columns. Per-caller QUAL via `INFO/VARQUAL[idx]`. Doubles as the catalog source for `COUNT_VARIANT`.
>>>>>>> remotes/origin/main

Exit criteria:
- Consensus sample-level VCF on the fixture carrying per-caller `FORMAT/{VAF,AD,DP,GT}` and cross-caller `INFO/{VAF_min,VAF_max,NCALLERS,DP_min,DP_median,DP_max,AD_REF_*,AD_ALT_*}`.
- Per-sample caller-comparison report (`pairwise_jaccard.tsv`, `af_pairwise_stats.tsv`, `summary.txt`) on the fixture.
- Per-caller `variants_count.h5ad` non-empty on the fixture for all enabled callers.
- `MIN_DP=200` and `MIN_BASE_QUALITY=30` removed from [bin/count_variant.py](bin/count_variant.py); `params.variant_min_dp` observable as a hard `INFO/DP` cutoff in the `BCFTOOLS_FILTER` output (records below threshold absent from the filtered VCF); `params.variant_min_vaf` observable as a row-count drop on per-cell rows whose `fraction_reads` is below the floor.

### Milestone 5 ; Reporting polish

Scope: surface the variant-calling layer in the metrics CSV, MultiQC report, and README. Active caller set at M5 design time is `{freebayes, bcftools, vardict, gatkhc}`; lofreq and deepvariant slot in transparently when re-enabled (the YAML and UpSet outputs scale by `params.variant_callers`).

Deliverables:

- [x] **COMB-468 — Identify and select metrics for MultiQC integration.** Captured in §4.9 (metric list, sources, descriptions) and §4.10 (UpSet plot). Decisions:
  - Thresholds blank in v1 (no truth set yet); revisit after benchmarking (§10.1) or first production rollout.
  - `report` column: per-caller + `consensus.*` → `aggregated`; `per_cell.<caller>.*` → `single`.
  - Consensus high-confidence cutoff: `INFO/NCALLERS >= 2`, **absolute** (not relative to active caller count); documented in the metric description.
  - UpSet plot in a sibling script (`bin/qc_variants.py`), not folded into `variants_to_metrics_parquet.py`.
- [x] **COMB-467 — Per-caller metrics rows in `metrics_reference.csv`.** 40 rows added in this branch (9 per-caller × 4 active callers + 4 `consensus.*`). All threshold columns blank per the decision above. Per-cell rows omitted (those metrics are deferred — see COMB-465 note).
- [x] **COMB-465 — `bin/variants_to_metrics_parquet.py` + `VARIANTS_TO_METRICS_YAML` module.** Per §4.9, staged in this branch. Module: [modules/variants_to_metrics_parquet/main.nf](modules/variants_to_metrics_parquet/main.nf). **Scope cut from spec draft:** the per-cell block (`per_cell.<caller>.*`) and the `consensus.ncallers_histogram` block are deferred to a follow-up; rationale in §4.9 "Deferred to a follow-up". Module input simplified from `(meta, parquet, consensus_vcf)` to `(meta, parquet)` — the parquet's `(chrom, pos, ref, alt)` × caller groups are equivalent to `INFO/NCALLERS` by §4.7 construction, so the VCF input would only re-derive the same signal.
- [x] **COMB-466 — MultiQC custom content.** Landed alongside `qc_variants.py` (2026-05-28). Two pieces:
  - **Per-caller + consensus metrics table.** Recon outcome: the Pair A "just works" claim was half-true — the variants modality auto-rendered, but as a single 40-column table. [bin/lib/metrics/report.py](bin/lib/metrics/report.py) `build_config` was extended to split the `variants` modality into two MultiQC sections by `metric_key.startswith("consensus.")`: `metrics_dna_variants_consensus` (4 headline metrics) and `metrics_dna_variants_per_caller` (36 per-caller metrics). Other modalities unchanged.
  - **UpSet PNG embed.** Mirror finding: `qc_amplicon`'s PNG was never actually mixed into MultiQC — it only lived on disk. Wired by emitting a `caller_upset_mqc.png` sibling from `QC_VARIANTS` (MultiQC auto-discovers `*_mqc.{png,jpg}` as standalone image sections) and mixing it into `ch_multiqc_files` at [main.nf](main.nf). The canonical `caller_upset.png` is still published under `samples/{sample}/qc/variants/` via the existing `variant_qc` channel.
- [x] **`bin/qc_variants.py` + `QC_VARIANTS` module.** Per §4.10. Sibling to `qc_amplicon.py`; reuses the `pip_pyyaml_upsetplot_numpy_pruned` container. [bin/lib/plotting/common.py](bin/lib/plotting/common.py) `plot_generic_upset_plot` gained three optional kwargs (backward compatible): `min_pct_show=None` to disable the subset-size floor, `highlight_min_degree` to highlight intersections of size ≥N via `UpSet.style_subsets(min_degree=...)` (replaces the qc_amplicon "all-flags-present" highlight), and `output_path` for an explicit filename. Module emits `caller_upset.png` (canonical published artifact), `caller_upset_mqc.png` (MultiQC sibling via `cp`), and `${sample}_variants_qc_metrics.yaml` (intersection counts as `upset.<sorted_callers>` keys). Wired in [subworkflows/dna.nf](subworkflows/dna.nf) off `PER_CELL_VARIANTS.out.variants_parquet`, the same channel feeding `VARIANTS_TO_METRICS_YAML`.
- [x] **COMB-469 — Integrate metrics extraction into `DNA`.** `VARIANTS_TO_METRICS_YAML.out.metrics_parquet` mixed into `DNA.out.metrics_parquet` at [subworkflows/dna.nf](subworkflows/dna.nf). `QC_VARIANTS.out.metrics_parquet` joined the same channel in the COMB-466 pass; `QC_VARIANTS.out.png` mixed into `variant_qc` (publishes to `samples/{sample}/qc/variants/`); new `variant_qc_mqc` emit carries the MultiQC sibling PNG into `ch_multiqc_files`.
- [ ] **Per-cell metrics follow-up ticket.** Not in the original ticket list; needs filing. Will add `per_cell.<caller>.{cells_with_calls,mean_variants_per_cell,median_variants_per_cell}` once the per-caller surface is validated in MultiQC. Requires joining `PER_CELL_VARIANTS.out.h5ad` (per-caller) back to a sample-level group before invoking a per-cell metrics script (or extending `variants_to_metrics_parquet.py` with additional h5ad inputs).
- [ ] **COMB-470 — README.** Germline-SNV scope, the new metrics + UpSet plot, removal of the legacy pseudo-bulk `bcftools` path (breaking change in release notes).

Acceptance criteria:
- `metrics_reference.csv` carries per-caller and `consensus.*` rows (and per-cell rows tagged `single`), thresholds blank.
<<<<<<< HEAD
- Aggregated MultiQC report shows **two** custom-content sections: `metrics_dna_variants_consensus` (high-confidence headline) and `metrics_dna_variants_per_caller` (per-caller `snv_count` / `singleton_count` / DP / VAF / het / hom). `consensus.ncallers_histogram` deferred — UpSet plot carries the same signal visually.
=======
- Aggregated MultiQC report shows **two** custom-content sections: `metrics_dna_variants_consensus` (high-confidence headline) and `metrics_dna_variants_per_caller` (per-caller `snv_count` / `pass_count` / `singleton_count` / DP / VAF / het / hom). `consensus.ncallers_histogram` deferred — UpSet plot carries the same signal visually.
>>>>>>> remotes/origin/main
- Aggregated MultiQC report carries one caller-agreement UpSet PNG per sample (via the `*_mqc.png` auto-discovery path). The canonical PNG is also published under `samples/{sample}/qc/variants/caller_upset.png`.
- README explicitly calls out germline/SNV scope and the legacy-path removal; release notes flag the latter as breaking.

---

## 10. Future Extensions

All extensions are additive; channel, layout, and subworkflow shape are already designed to accommodate them.

### 10.1 Benchmarking

Concordance against a user-supplied truth VCF, using RTG Tools `rtg vcfeval` (GA4GH small-variant comparison; Apache 2.0). hap.py was the historical alternative but has had no substantive release since ~2019; `rtg vcfeval` is the maintained equivalent.

- New module `modules/benchmark_variants.nf` with `PREPARE_REFERENCE_SDF`, `BENCHMARK_VARIANTS`, `MERGE_BENCHMARK_PARQUET`; helper `bin/benchmark_summary.py`.
- New params: `truth_vcf`, `truth_vcf_tbi`, `confident_regions_bed` (defaults to `amplicon_bed`), `benchmark_tool` (`vcfeval` \| `isec`).
- New samplesheet column: `truth_vcf` (per-sample override).
- New container: `rtg-tools=3.12.1` (Wave build).
- Per `(sample, caller)` outputs: `summary.txt`, `tp/fp/fn.vcf.gz`, `weighted_roc.tsv.gz`; merged into `experiment/variant_benchmark.parquet` with columns `tp, fp, fn, precision, recall, f1, truth_total, query_total`.
- Surface per-caller precision/recall/F1 in `variants_to_metrics_parquet.py` and ROC TSVs in MultiQC custom content.

### 10.2 Indel support
- Drop `--no-indels` / `--skip-indels` / equivalent from each caller.
- Add LoFreq `indelqual --dindel` pre-step (writes to a work-dir BAM copy ; stages are read-only).
- Drop `-v snps` from `BCFTOOLS_NORM`.
- Add `type is_a {SNV, INS, DEL, MNV}` column to `VCF_TO_PARQUET` (`VCF_TO_COUNTS_CSV` no longer exists — see §4.6).
- Introduce `bin/count_variant_indel.py` with an INS/DEL branch using `pysam.PileupRead.indel` / `is_del`; feature key becomes `CHROM:POS:REF>ALT`.

### 10.3 Somatic support
- Add `params.variant_mode is_a {germline, somatic}` (default `germline`) + optional samplesheet override column.
- Per-caller somatic args (via `task.ext.args_somatic` in `conf/modules.config`):
  - freebayes: `--pooled-continuous --min-alternate-fraction 0.01 --min-alternate-count 2 --min-coverage 50 --haplotype-length 0`
  - LoFreq: add `--min-alt-bq 20 --sig 0.01`
  - bcftools: `--ploidy 1` (pseudo-haploid low-VAF trick)
  - VarDict: `-f 0.01`; consider `var2vcf_paired.pl` if tumor/normal pairing arrives
- Document the caveat: "somatic" on pseudo-bulk amplicon ≠ tumor-normal calling.
`TODO: add GATK Mutect to the somatic callers`

### 10.4 Error correction

- For error correction two tools are considered for this workflow: UMI-tools and fgbio.
- Investigate how these two tools can be used to collapse reads and error correction.
- Error correction for the cell barcode based on base quality; A probablistic model; UMI-tools also supports different algorithm in grouping metrics

---

## 11. Known Risks

| Risk | Mitigation |
|---|---|
| Removing the legacy `PILEUP_READS` path drops an output some downstream processes still reads |  Flag anything found in the Milestone 1 PR. |
| VarDict produces zero calls silently on bad/empty BED | A post-call record count > 0 on the fixture or maybe a simple `wc -l` for rows not starting `#` (i.e. VCF header) |
| DeepVariant container is ~7 GB and slow to pull on a cold cache or first run | Pre-pull on CI workers; According to documents, we should pin a specific version `google/deepvariant:1.6.1` and not `:latest`. |
| DeepVariant runtime can vary and it can dominate pipeline run time on large panels | Create a `num_shards` via `task.cpus`; document that DeepVariant can be dropped from `params.variant_callers` for rapid-iteration runs. |
| GATK HaplotypeCaller can run out of memory inside the container due to JVM default values. | Always set `--java-options "-Xmx${task.memory.toGiga()}g"` in the process body |
| Upstream tools disappearance (especially LoFreq, see §7.2) | Per-caller decoupling via `params.variant_callers`; individual callers can be dropped without an architectural change. Re-review §7.2 each major release. |
| Container availability for the M2/M3 callers | v1 sources callers from BioContainers pinned-hash images (§7). If a tool lacks a usable BioContainers tag, fall back to a Wave build for that tool only. |
| A developer introduced error/bug that a Channel joins on `meta` instead of `(meta, caller)` | Convention: every downstream join on variant-calling channels uses `by: [0, 1]`. Add a comment in the subworkflow header. `TODO: is it possible to generate a verification step for this to automate this potential bug?`|
| Normalization left-alignment surprises for SNVs or even InDels | Document Document Document. Raw VCFs remain available via `raw_vcf`. |
| User sets `variant_callers: []` expecting it to disable variant calling | `nextflow_schema.json` enforces `minItems: 1`; validation error message points to `run_variant_calling` as the disable switch. |
| Per-sample `run_variant_calling: true` set on a non-DNA row | Ignored silently in the filter (since this logic will only affect the `DNA_PROCESSING`); `validateParameters()` does not need to enforce this. `TODO: Document in the samplesheet schema description.` |

---

## 12. Open Decisions

1. `.vcf.gz` vs. `.bcf`? *(proposed: `.vcf.gz`)*
2. Trust caller-native `INFO/AF`, or always recompute from `AD` in parquet? *(proposed: prefer native, compute when absent)*
3. Where does `amplicon_bed` come from? An existing asset accompanying `FASTA`, a new required input, or derived from `amplicon_fasta` headers? *(proposed: derive as fallback from FASTA header; `params.amplicon_bed` overrides)*
4. QUAL cutoff defaults per caller ; accept the values in §6 or tune empirically once a truth set is available?
5. `params.run_variant_calling` default ; `true` (current proposal, opt-out) or `false` (opt-in, since variant calling is the heaviest part of `DNA_PROCESSING`)?
6. Should a per-sample `run_variant_calling: true` on a non-DNA (RNA) row warn at parse time, or stay silent?

---

## Appendix A ; Files that will be affected (modify, create, delete)

### To modify
- [subworkflows/dna_processing.nf](subworkflows/dna_processing.nf) ; calls `VARIANT_CALLING`; produces `variant_vcf`, `variant_stats`, `variant_parquet`. Remove `PILEUP_READS`!
- [modules/dna_processing.nf](modules/dna_processing.nf) ; **delete the `PILEUP_READS` process** ([modules/dna_processing.nf:77-125](modules/dna_processing.nf#L77-L125)).
- [main.nf](main.nf) ; new entries under `samples/{sample}/DNA/variants/{caller}/`; mix caller stats into `ch_multiqc_files`. Remove any code pointing at legacy `snps.csv` / pseudo-bulk VCF outputs.
- [nextflow_schema.json](nextflow_schema.json) ; new params per §6.1, including `run_variant_calling` (boolean) and `variant_callers` with `minItems: 1`. Remove any legacy variant-calling schema entries tied to the removed path.
- [assets/schema_input.json](assets/schema_input.json) ; new optional column (`run_variant_calling`).
- [assets/metrics/metrics_reference.csv](assets/metrics/metrics_reference.csv) ; per-caller metric rows; delete rows describing legacy pseudo-bulk metrics if any exist.
- `conf/modules.config` ; `task.ext.args` for each caller.
- [README](README.md) ; remove the "quite poor pseudo-bulk bcftools" disclaimer; replace with the six-caller description and scope notes; add a concise summary of the variant calling workflow.

### To create
- `subworkflows/variant_calling.nf`
- `modules/samtools/faidx/main.nf` (`SAMTOOLS_FAIDX`); `modules/gatk4/createsequencedictionary/main.nf` (`GATK4_CREATESEQUENCEDICTIONARY`, M2); `modules/gatk4/updatevcfsequencedictionary/main.nf` (`GATK4_UPDATEVCFSEQUENCEDICTIONARY`, M2 — lofreq reheader)
- `modules/bcftools/norm/main.nf` (`BCFTOOLS_NORM`); `modules/bcftools/stats/main.nf` (`BCFTOOLS_STATS`)
- `modules/tabix/main.nf` (`TABIX`)
- `modules/freebayes/main.nf` (`FREEBAYES`); `modules/lofreq/call/main.nf` (`LOFREQ_CALL`, M2); per-tool `modules/<caller>/main.nf` (or per-subcommand `modules/<caller>/<op>/main.nf` for multi-op tools) for the remaining M2/M3 callers
- `bin/vcf_to_parquet.py`
- `bin/variants_to_metrics_parquet.py`

### To delete
- Purge any helper script or fixture tied exclusively to the legacy `PILEUP_READS` → `snps.csv` path.

### Not needed in v1 (reserved for future)
- `modules/benchmark_variants.nf` and `bin/benchmark_summary.py` ; only required by the benchmarking extension (§10.1).
- `bin/count_variant_indel.py` ; only required by the indel extension.
- `bin/vcf_to_counts_csv.py` ; removed in M4 ; `bin/vcf_to_parquet.py` is the only catalog-extraction script (§4.7).

## Appendix B ; Other existing workflows for variant callers

The list of existing variant calling workflows is extensive. And the following list will grow. The idea is to understand where we can use the existing knowledge in the OpenSource community when it comes to the difficult task of including variant callers. Mainly to address the risks in describe in §11 and §7

### [nf-core/sarek v3.4.0](https://nf-co.re/sarek/3.4.0)

> A Nextflow pipeline workflow designed to detect variants on whole genome or targeted sequencing data. Initially designed for Human, and Mouse, it can work on any species with a reference genome. Sarek can also handle tumour / normal pairs and could include additional relapses.

### Lessons Learned:

The choice of variant caller from Sarek, especially in Germline, can help us choose our tools. As we can see from the risk categories in §11 and §7, we have an inherent risk of drowning ourselves in a long list of variant callers that either poorly maintined or not maintained at all.

The branching logic in Uniflow is inspired by Sarek, since it is dealing with both Somatic and Germline.

### [nf-core/raredisease v1.1.1](https://nf-co.re/raredisease/1.1.1)

> nf-core/raredisease is a best-practice bioinformatic pipeline for calling and scoring variants from WGS/WES data from rare disease patients. This pipeline is heavily inspired by [Clinical Genomics - MIP](https://github.com/Clinical-Genomics/MIP).

### Lessons Learned:

The lean structure, and the choice of limited variant callers for germline rare variant calling is obvious. Ultimately we should aim for including only a few variant callers that are industry standards.

The choice of DeepVariant was inspired by raredisease workflow. And the options for pileup from MIP.

### [BALSAMIC](https://github.com/Clinical-Genomics/BALSAMIC)

> A somatic only variant caller with options to use Germline calls for filtering. It is an ISO accredited workflow currently in use in Karolinska Hospital (Disclaimer: Hassan Foroughi was the lead developer of this workflow)

### Lessons Learned:

A somatic only workflow which heavily relies Sentieon, a commercially available product, although the opensource variant callers such as Mutect2 and VarDict also exist as a choice.

---

## Appendix C ; Cell barcode collapsing for error correction

This section describes the left out strategy of error correction. There are two main tools we can use here. This appendix covers the main points and workflow.

### fgbio

`fgbio` is a collection of tools for grouping, collapsing, and consensus calling of the reads based on UMI and/or cell barcodes. Their current best practice is a simple [workflow](https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md):

```
fastqToBam to 
minimap2 to # First mapped BAM
GroupReadsByUMI to 
CallMolecularConsensus to 
minimap2 # Second mapped BAM
```

The two BAM files generated can be used for bulk variant callers to call SNV and InDels.

### UMI-tools

TBD