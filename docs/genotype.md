<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Reports](reports.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Understand genotype outputs

Uniflow pipeline discovers SNVs and indels in the DNA amplicon library and genotypes individual cells at each discovered variant. This page describes where the genotypes are, how they are called, and what each field means.

> [!WARNING]
> Variant-calling functionality is under active development and will continue to evolve. The current release supports SNV and indel detection, together with per-cell genotype and likelihood assignments. While these capabilities have undergone initial testing, validation across the full range of samples, experimental conditions, and use cases is ongoing. Variant-calling results should therefore be interpreted with appropriate scientific judgment.

## Where to find genotypes

Variant outputs are produced when the selected product enables variant calling (see [Select a product](inputs.md#select-a-product)).

| Output | Contents |
| --- | --- |
| `samples/<sample_name>/counts/counts.h5mu` | Modalities `snv` and `indel`: cell-by-position and cell-by-indel matrices with genotype (`GT`), genotype quality (`GQ`) and read-count layers. `counts.zarr/` holds the same data. |
| `samples/<sample_name>/dev/variant/variants_snv.parquet` | One row per cell, position and observed allele, with read counts, qualities, genotype likelihoods and caller support. |
| `samples/<sample_name>/dev/variant/variants_indel.parquet` | The same columns for insertions and deletions. |
| `samples/<sample_name>/DNA/<sample_name>.consensus.vcf.gz` | Variants discovered from the pooled reads of all cells, with one sample column per variant caller and a `.tbi` index. Contains no per-cell data. |
| `samples/<sample_name>/metrics/metrics.csv` | Discovery and attribution metrics, in the rows whose `modality` is `variants`. |

Start with `counts.h5mu` for downstream analysis. Use the parquet tables when you need allele-level detail or genotype likelihoods. As `dev/` outputs, their columns may change between releases.

## Before you start

- **`GT` -1 means not genotyped, not reference.** Homozygous reference is `GT` 1. Low-quality genotypes keep their `GT`, so also filter on `GQ`. See the table below in section [Count Matrices](#count-matrices) for further explanation.
- **Select called cells first.** All modalities in `counts.h5mu` share one cell axis that holds every barcode seen in any modality. Select cells with, for example, `mdata.obs["amplicon:is_cell"]`.
- **Positions without a column are unknown.** The `snv` modality has a column only for discovered positions that at least one cell read reached. A position without a column was not examined, not found to be reference.
- **Indel non-carriers are genotyped.** A cell that covers an indel without a read carrying it is genotyped from its reference reads, normally as homozygous reference (`GT` 1).
- **Only discovered variants are genotyped, from downsampled reads.** A variant with too little support in the pooled reads has no per-cell genotype, and read depths are capped at 100 reads per cell and amplicon by default.

## How genotypes are called

1. **Alignment:** Reads are aligned to the amplicon FASTA. This full-depth alignment is the BAM in `DNA/`.
2. **On-target reads:** A read is on-target when it is a primary alignment with a mapping quality of at least 20 that starts at the most common read start or ends at the most common read end of its amplicon. On-target reads are counted per amplicon to call cells.
3. **Downsampling:** By default, variant analysis uses at most 100 on-target reads per cell and amplicon, sampled with a fixed seed from called cells only. The amplicon counts and the BAM in `DNA/` keep full depth, so read depths in the genotype outputs are lower by design. With downsampling disabled, steps 4 and 5 use the full-depth alignment.
4. **Discovery:** FreeBayes, bcftools, GATK HaplotypeCaller and VarDict call variants on the pooled reads of all cells. Their calls are split into one ALT allele per record, left-aligned, and kept when at least 20 reads support the ALT allele at a VAF of at least 0.0001 (0.001 for FreeBayes and VarDict). The calls of all callers are merged into the consensus VCF; a variant reported by any one caller is genotyped.
5. **Per-cell counting:** At every discovered position, Uniflow counts each cell's reads per allele. Only on-target reads are counted, whether or not downsampling is enabled. Positions covered by fewer than 11 reads across all cells are skipped.
6. **Genotyping:** Each cell is genotyped from its reference and ALT read counts, as described below.
7. **Labelling:** Weak observations and low-quality genotypes are labelled, never removed.

SNVs, insertions and deletions are genotyped. Multi-nucleotide variants (MNVs) appear in the consensus VCF but do not receive a per-cell genotype in the current pipeline implementation.

### Genotype model

For each cell, Uniflow compares three diploid genotypes (homozygous ref, homozygous alt, heterozygous) using the reads supporting the reference (`ref_reads`) and the ALT allele (`alt_reads`). Reads showing any other allele, or no base call (`N`), are not used.

| Genotype | `gt` | Expected ALT read fraction |
| --- | --- | --- |
| Homozygous reference | `0/0` | 0.00045, the assay error rate. For SNVs, raised when the ALT reads have low base quality. |
| Heterozygous | `0/1` | 0.5 |
| Homozygous ALT | `1/1` | 0.99955 |

Read counts follow a beta-binomial distribution with an overdispersion of 0.091, calibrated once on a reference sample with known genotypes. `gq` can be optimistic for samples that are more overdispersed than that reference.

- `gt` is the most likely genotype.
- `pl_homref`, `pl_het` and `pl_homalt` are Phred-scaled genotype likelihoods, shifted so that the called genotype is 0.
- `gq` is the second-smallest PL, capped at 99: how much more likely the called genotype is than the next best.
- A genotype with `gq` below 30 is labelled `lowGQ`. In `counts.h5mu`, it keeps its `GT` code; only `GQ` marks it.
- `outlier_lod` is the log10 likelihood ratio of an unconstrained ALT read fraction against the called genotype. Positive values mean the reads fit none of the three genotypes well, as expected from doublets, ambient DNA or unequal allele amplification. It does not change `gt` or `gq`.

Read depth limits `gq`. A homozygous call with no reads of the other allele reaches GQ 30 at about 21 reads. A balanced heterozygous call reaches GQ 30 at 4 reads but levels off near GQ 50 however deep the cell is sequenced. A single Q40 ALT read among 100 reads lowers a homozygous-reference call to GQ 29, so `lowGQ` `0/0` genotypes are either shallow or carry a few ALT reads, so we treat them as undetermined.

### ALT allele of a genotype

- **SNVs:** Each cell gets one genotype per position. Its ALT allele is the non-reference base with the most reads in that cell; ties go to the alphabetically first base. At a position with several ALT bases, cells can be genotyped against different bases, and **a cell's ALT is not necessarily an allele a caller reported**. Check `alt` and `callers` on the genotyped row.
- **Indels:** Each cell with a read spanning the anchor base of a called indel receives one genotype for that indel. A cell without a read carrying the indel is genotyped from its reference reads.

## Per-cell tables

`variants_snv.parquet` and `variants_indel.parquet` share one schema and can be concatenated.

### Rows

In `variants_snv.parquet`, each cell has, for each discovered position that it covers:

- one `REF` row counting its reference reads, if it has any
- one `ALT` row for each non-reference base it shows, including `N`

The cell's genotype is on exactly one of these rows, marked by `genotyped`: the `ALT` row of the base the cell was genotyped against, or the `REF` row when the cell shows no non-reference base other than `N`. The genotype columns (`variant_name`, `gt_filter`, `gt`, `gq`, the PLs, `outlier_lod` and `zygosity`) are empty on all other rows.

In `variants_indel.parquet`, each row is one cell and one called indel, for every cell with a read spanning the indel's anchor base; `alt_reads` is 0 for cells without a read carrying it. A called indel that no cell spans has a single row with an empty `barcode`.

To get one genotype per cell and site, keep the rows where `genotyped` is true. To keep confident genotypes only, also require `gt_filter` to be `PASS`:

```python
import pandas as pd

table = pd.read_parquet("samples/<sample_name>/dev/variant/variants_snv.parquet")
genotypes = table[table["genotyped"] & (table["gt_filter"] == "PASS")]
```

or using Polars:
```python
import polars as pl

genotypes = (
    pl.scan_parquet("samples/<sample_name>/dev/variant/variants_snv.parquet")
    .filter(pl.col("genotyped") & (pl.col("gt_filter") == "PASS"))
    .collect()
)
```

### Columns

| Column | Description |
| --- | --- |
| `sample` | Sample name. |
| `barcode` | Cell barcode. |
| `kind` | `SNV` or `INDEL`. |
| `record_type` | `REF` for reference-read rows (SNVs only), otherwise `ALT`. |
| `feature` | Amplicon ID from the amplicon FASTA. |
| `pos_local_0based` | 0-based position in the amplicon. For indels, the position of the reference base before the event. |
| `chrom`, `pos` | Experimental chromosome and 1-based genome position. Empty unless the amplicon FASTA carries genome coordinates (see [Coordinates](#coordinates)). |
| `ref` | SNV: reference base. Indel: reference allele including the base before the event, as in a VCF. |
| `alt` | SNV: observed base, which may be `N`. Indel: ALT allele as in a VCF. Empty on `REF` rows. |
| `variant_name` | Genotype label: `homo<ref>><alt>` for homozygous ALT (for example `homoC>T`), `homo<ref>><ref>` for homozygous reference, `hete` for heterozygous, and `lowqual` for `lowGQ` genotypes. |
| `total_reads` | SNV: reads with a base call at the position, including `N`. Indel: all reads spanning the base before the event. |
| `ref_reads` | Reads supporting the reference at the site in this cell. For indels, reads with no indel at the site. |
| `alt_reads` | Reads supporting `alt`. Empty on `REF` rows. |
| `vaf` | `alt_reads / total_reads`. Empty on `REF` rows. |
| `mean_bq` | Mean base quality of the row's reads, averaged as error probabilities and reported on the Phred scale. Empty for indels. |
| `mean_mq` | Mean mapping quality of the row's reads. |
| `read_filter` | `PASS`, or the reasons the row's reads are weak, joined by `;`: `lowDP` (`total_reads` below 5), `lowVAF` (the row's read fraction below 0.01), `lowBQ` (`mean_bq` below 15) and `lowMQ` (`mean_mq` below 30). Does not affect the genotype. |
| `gt_filter` | `PASS` when `gq` is at least 30, otherwise `lowGQ`. |
| `is_no_call` | SNV: true when `alt` is `N`. Empty for indels. |
| `cells_supporting_alt` | Cells with at least 5 reads at the site and a read fraction of at least 0.01 for this allele. Empty on `REF` and `N` rows. |
| `cells_total_at_site` | Cells with at least 5 reads at the site. |
| `gt` | Most likely genotype: `0/0`, `0/1` or `1/1`. Present for `lowGQ` genotypes too. |
| `gq` | Genotype quality, from 0 to 99. |
| `pl_homref`, `pl_het`, `pl_homalt` | Phred-scaled genotype likelihoods; the called genotype is 0. |
| `outlier_lod` | Fit of the reads to any diploid genotype; see [Genotype model](#genotype-model). |
| `zygosity` | `homref`, `het` or `hom` (homozygous ALT). Empty for `lowGQ` genotypes. |
| `callers` | Variant callers that reported this exact allele. Empty for alleles no caller reported. |
| `n_callers` | Number of callers in `callers`. |
| `genotyped` | True on the row carrying the cell's genotype for the site. |

## Count matrices

Load the counts, select called cells and count the carrier cells at each measured SNV position:

```python
import mudata as mu
import numpy as np

mdata = mu.read_h5mu("samples/<sample_name>/counts/counts.h5mu")
cells = mdata.obs_names[mdata.obs["amplicon:is_cell"]]

snv = mdata["snv"][cells]
gt, gq = snv.layers["GT"], snv.layers["GQ"]
carrier_cells = (np.isin(gt, [2, 3]) & (gq >= 30)).sum(axis=0)
```

Both modalities encode genotypes the same way:

| `GT` | Meaning |
| --- | --- |
| `-1` | Not genotyped: the cell does not cover the position, or, for indels, no read spanning the anchor base shows the reference or this indel (`REF` + `AD` is 0) |
| `1` | Homozygous reference |
| `2` | Heterozygous |
| `3` | Homozygous ALT |

`0` is never written. Every genotype keeps its `GT` code whatever its `GQ`, so select confident genotypes with `GT` in 1, 2 or 3 and `GQ` of at least 30. `GQ` is `NaN` exactly where `GT` is -1; a `GQ` of 0 is a real value.

Short definitions of `X`, depth and `GT` are also stored in each object, in `mdata["snv"].uns["snv"]` and `mdata["indel"].uns["indel"]`. In `mdata.obs`, the per-cell columns below carry the `snv:` or `indel:` prefix.

### `snv`: cells by amplicon positions

Columns are the discovered positions that at least one cell read reached, named `<amplicon_id>_<amplicon_pos>`. Other amplicon positions have no column. `uns["declared_target_sites"]` lists the positions declared as targets in the amplicon FASTA and whether each was measured.

`X` is the non-reference read fraction: `(DP - reads matching ref_base) / DP`, counting every non-reference base.

| Layer | Contents |
| --- | --- |
| `AD_A`, `AD_C`, `AD_G`, `AD_T` | Reads showing each base |
| `AD_N` | Reads showing `N` or another base outside A, C, G and T |
| `DP` | `AD_A + AD_C + AD_G + AD_T`. Excludes `AD_N`, so `DP + AD_N` equals `total_reads` in the parquet table. |
| `GT`, `GQ` | Genotype code and genotype quality |

`GT` does not identify the ALT base. At positions with several ALT bases, use the `AD_*` layers or the `alt` column of the parquet table.

| `var` column | Description |
| --- | --- |
| `amplicon_id` | Amplicon ID |
| `amplicon_pos` | 1-based position in the amplicon |
| `ref_base` | Reference base |
| `is_measured` | Always true, since only measured positions have columns. A discovered position with too few reads has no column. |
| `is_target_site` | Experimental and reserved for position matching `snp_pos` in the amplicon FASTA header |
| `is_candidate_variant` | A caller reported a variant at this position |
| `n_callers` | Callers reporting any allele at this position. In the parquet tables, `n_callers` counts callers per allele instead. |
| `pseudobulk_depth` | Median depth reported by the callers |
| `pseudobulk_nonref_count` | Median ALT reads reported by the callers, summed over the position's ALT alleles |
| `pseudobulk_nonref_vaf` | `pseudobulk_nonref_count / pseudobulk_depth`. Can exceed 1 at positions with several ALT alleles. |
| `pseudobulk_major_nonref` | ALT allele with the most pooled ALT reads |
| `n_cells_covered` | Cells with `DP` above 0 |
| `n_cells_with_nonref_support` | Cells with at least one non-reference read |

| `obs` column | Description |
| --- | --- |
| `total_amplicon_depth` | Sum of `DP` over the cell's measured positions. A read spanning several positions is counted at each. |
| `n_sites_covered` | Positions with `DP` above 0 |
| `mean_site_depth`, `median_site_depth` | Mean and median `DP` over covered positions |
| `n_sites_with_nonref_support` | Covered positions with at least one non-reference read |
| `fraction_sites_with_nonref_support` | `n_sites_with_nonref_support / n_sites_covered` |

### `indel`: cells by called indels

Columns are the called insertions and deletions, one per allele. An indel is anchored at the reference base before the event, as in a VCF.

`X` is `AD / DP`.

| Layer | Contents |
| --- | --- |
| `AD` | Reads carrying this indel |
| `REF` | Reads with no indel at the anchor base |
| `OTHER_INDEL` | Reads with a different indel at the anchor base, or whose anchor base is deleted by a larger event |
| `DP` | All reads spanning the anchor base: `AD + REF + OTHER_INDEL` |
| `GT`, `GQ` | Genotype code and genotype quality |

To count the cells covering an indel without carrying it, use `DP` above 0 with `AD` equal to 0, not `GT`. `DP` includes reads that do not span a whole deletion, so indel VAF is conservative.

| `var` column | Description |
| --- | --- |
| `event_id` | `<amplicon_id>:<event_type>:<start_pos>:<end_pos>:<alt_seq>`. Column names equal `event_id`, with `@<pos_local_0based>:<ref>><alt>` appended when two alleles share one. |
| `amplicon_id` | Amplicon ID |
| `event_type` | `ins` or `del` |
| `start_pos`, `end_pos` | 1-based amplicon positions. Deletion: first and last deleted base. Insertion: the two bases the insertion lies between. |
| `ref_seq`, `alt_seq` | Deleted or inserted bases without the anchor base; `-` when empty |
| `event_len` | Length: positive for insertions, negative for deletions |
| `pos_local_0based`, `ref`, `alt` | The indel as in a VCF: 0-based anchor position and alleles including the anchor base. These match `feature`, `pos_local_0based`, `ref` and `alt` in `variants_indel.parquet`. |
| `n_callers` | Callers reporting this indel |
| `is_candidate_indel` | Always true |
| `pseudobulk_depth`, `pseudobulk_nonref_count`, `pseudobulk_nonref_vaf` | Median depth and median ALT reads reported by the callers, and their ratio |
| `indel_feature` | `<amplicon_id>:<pos_local_0based>:<ref>><alt>` |
| `n_cells_covered` | Cells with `DP` above 0 |
| `n_cells_with_indel_support` | Cells with `AD` above 0 |
| `event_qc_pass` | True when at least one cell carries the indel, false when cells cover the anchor base but none carries it, empty when no cell covers it |
| `is_measured` | At least one cell covers the anchor base (`n_cells_covered` above 0) |
| `is_target_site` | Always empty for indels |

| `obs` column | Description |
| --- | --- |
| `total_amplicon_depth` | Sum of `DP` over the distinct anchor bases the cell covers |
| `n_indel_events_covered` | Indels with `DP` above 0 |
| `n_indel_events_with_support` | Indels with `AD` above 0 |
| `total_indel_support` | Sum of `AD` |
| `mean_event_depth`, `median_event_depth` | Mean and median `DP` over covered indels |
| `fraction_indel_events_with_support` | `n_indel_events_with_support / n_indel_events_covered` |

## Consensus VCF

`<sample_name>.consensus.vcf.gz` lists the discovered variants. It describes the pooled reads of all cells.

- `#CHROM` is the amplicon ID and `POS` is the 1-based position in the amplicon.
- Each record has one ALT allele; multi-allelic sites are split into several records.
- Each variant caller has its own sample column, named after the caller. Column order can differ between runs, so select callers by name. A caller that did not report the allele has missing values.

| FORMAT field | Description |
| --- | --- |
| `GT` | The caller's genotype for the pooled reads. Not a cell genotype. |
| `DP` | Read depth reported by the caller |
| `AD` | Reads supporting the reference and ALT alleles |
| `VAF` | ALT reads divided by the sum of `AD` |
| `FILTER` | The caller's own filter status |
| `VARQUAL` | The caller's variant quality |

| INFO field | Description |
| --- | --- |
| `NCALLERS` | Callers that reported the allele |
| `DP_min`, `DP_median`, `DP_max` | Depth across callers |
| `AD_REF_min`, `AD_REF_median`, `AD_REF_max` | Reference reads across callers |
| `AD_ALT_min`, `AD_ALT_median`, `AD_ALT_max` | ALT reads across callers |
| `VAF_min`, `VAF_max` | ALT read fraction across callers |

The `QUAL` and `FILTER` columns are carried over when the callers' records are merged; they are not a consensus score. Use `NCALLERS` and each caller's `VARQUAL` and `FILTER` instead.

## Coordinates

This section is experimental FASTA format fields. It comes with no support for this beta release. These are reserved fields in the FASTA reference file that we are planning to utilize in future releases, and otherwise should be ingored.

| Output | Field | Coordinate |
| --- | --- | --- |
| Parquet tables | `pos_local_0based` | 0-based, amplicon |
| Parquet tables | `pos` | 1-based, genome |
| `snv` modality | `amplicon_pos` and column names | 1-based, amplicon |
| `indel` modality | `start_pos`, `end_pos` and column names | 1-based, amplicon |
| `indel` modality | `pos_local_0based` | 0-based, amplicon |
| Consensus VCF | `POS` | 1-based, amplicon |

For an SNV, `amplicon_pos` and the VCF `POS` both equal `pos_local_0based + 1`.

Genome coordinates are filled in only when each amplicon FASTA header gives the chromosome (`chrom`) and the genome position of the first amplicon base (`g0`), with the amplicon sequence on the forward strand:

```text
>AMPLICON_1 chrom=chr11 g0=67734055 snp_pos=67734444
```

`pos` is then `g0 + pos_local_0based`. The optional `snp_pos` marks the designed target position (`is_target_site`). Write the keys in lowercase and give a single `snp_pos`; other forms pass input validation but are ignored.

## Variant metrics

Rows whose `modality` is `variants` in `samples/<sample_name>/metrics/metrics.csv`:

| Metric | Description |
| --- | --- |
| `snvs_called`, `indels_called`, `mnvs_called` | Distinct ALT alleles in the consensus VCF |
| `snvs_attributed`, `indels_attributed`, `mnvs_attributed` | Called alleles with at least one supporting read in a cell. Called but unattributed variants have no genotypes; MNVs are never attributed. |
| `consensus.total_sites` | Distinct SNV alleles in the consensus VCF |
| `consensus.high_confidence_count`, `consensus.singleton_count`, `consensus.fraction_high_confidence` | SNV alleles reported by at least 2 callers, by exactly 1 caller, and the fraction reported by at least 2 |
| `consensus.cell_support_sites` | Non-reference alleles at discovered SNV positions, whether a caller reported them or not, with a non-empty `cells_supporting_alt` |
| `consensus.median_cells_supporting_alt`, `consensus.multicell_support_count`, `consensus.singlecell_support_count`, `consensus.fraction_multicell_support` | Distribution of `cells_supporting_alt` over those alleles: its median, the alleles supported by at least 2 cells or by exactly 1, and the fraction supported by at least 2 |
| `<caller>.snv_count`, `<caller>.singleton_count` | SNV alleles reported by the caller, and by that caller only |
| `<caller>.mean_dp`, `<caller>.median_dp`, `<caller>.mean_vaf`, `<caller>.median_vaf` | Pooled depth and VAF of the caller's SNVs |
| `<caller>.het_count`, `<caller>.hom_count` | The caller's SNVs with a pooled VAF above 0.05 and below 0.8, or of at least 0.8. These are VAF bands, not genotypes. |
| `upset.<callers>`, `upset_indel.<callers>` | SNV and indel alleles reported by exactly that combination of callers |

## Settings

| Parameter | Default | Effect |
| --- | --- | --- |
| `--product_id` | Required | Products ending in `-no-variant` skip variant calling. |
| `--variant_callers` | `freebayes,bcftools,gatkhc,vardict` | Callers used for discovery. Any subset may be given. |
| `--per_cell_per_amplicon_downsample` | `true` | Caps the reads per cell and amplicon used for variant analysis. With `false`, variant analysis uses all reads at full depth, including reads from barcodes not called as cells. |
| `--downsample_reads_per_cell_amplicon` | `100` | Maximum reads per cell and amplicon |
| `--downsample_seed` | `42` | Seed for the read sampling |

The other thresholds on this page, including the discovery floors, the genotype model and the GQ threshold, are set in `conf/variant.config`.

## Documentation sections

1. [Getting started](docs/getting-started.md)
2. [Installation](docs/installation.md)
3. [Prepare inputs](docs/inputs.md)
4. [Run Uniflow](docs/running.md)
5. [Understand outputs](docs/outputs.md)
6. [Output definitions](docs/specs/outputs.md)
7. [QC report files](docs/reports.md)
8. [Genotyping](docs/genotype.md)
9. [Troubleshooting and support](docs/support.md)