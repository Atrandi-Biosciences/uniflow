# uniflow output overview

This document describes the high-level structure of a `uniflow` output folder for the single-cell RNA + DNA assay. It explains what each major output contains, what the files represent, and which files are most useful for reviewing data quality or continuing downstream analysis.

## Top-level structure

```text
uniflow_v0.5.0/
├── experiment/
├── libraries/
├── metadata/
└── samples/
```

| Path | Description |
| --- | --- |
| `experiment/` | Run-level summary outputs across all libraries and samples. Start here for an overview of the full run. |
| `libraries/` | Per-library outputs. A library usually corresponds to one sequencing library or lane-specific FASTQ pair, such as `lib1_RNA_L001` or `lib2_DNA_L002`. |
| `metadata/` | Input sample metadata and expanded demultiplexing information used by the pipeline. |
| `samples/` | Per-sample outputs after RNA and DNA libraries have been assigned to biological samples and processed together. |

## Recommended files to review first

| File | Why it matters |
| --- | --- |
| `experiment/experiment.html` | Main run-level HTML report. Use this as the starting point for reviewing the experiment. |
| `experiment/metrics.csv` | Combined metrics table across libraries and samples. Useful for comparing all outputs in one place. |
| `samples/<sample_id>/metrics/metrics.csv` | Per-sample metrics for RNA gene expression and DNA amplicon performance. |
| `samples/<sample_id>/counts/counts.h5mu` | Main combined count object for downstream single-cell analysis. |
| `samples/<sample_id>/qc/` | Per-sample QC plots for gene expression, amplicon, and carryover behavior. |
| `libraries/<library_id>/qc/barcode/` | Per-library barcode QC plots. Useful for diagnosing barcode parsing and library-level quality. |

## `experiment/`

```text
experiment/
├── experiment.html
├── metrics.csv
└── qc_support.tar.gz
```

### `experiment.html`

The top-level HTML report summarizes the full pipeline run. It is intended as the primary human-readable report for customers and support teams.

Use it to review:

- Overall run quality.
- Library-level barcode performance.
- Sample-level RNA and DNA metrics.
- QC plots and summary tables collected across the experiment.

### `metrics.csv`

This is a long-format metrics table combining library-level and sample-level metrics.

Columns include:

| Column | Meaning |
| --- | --- |
| `source_id` | The library or sample the metric belongs to. Examples: `lib1_RNA_L001`, `sample1`. |
| `library_type` | The library type, usually `RNA` or `DNA`. |
| `modality` | The analysis module that produced the metric, such as `barcode`, `gene_expression`, or `amplicon`. |
| `metric_key` | Name of the metric. |
| `value` | Metric value. |

This file is useful for programmatic review, custom dashboards, and comparing samples or libraries across a run.

### `qc_support.tar.gz`

This file contains any relevant data to be shared with our customer support if you need help with troubleshooting your run.

The files included **DO NOT** provide sensitive information regarding the underlying data but only QC relevant information.

## `metadata/`

```text
metadata/
├── sample_sheet.csv
```

### `sample_sheet.csv`

This file records the input FASTQ files and their sample assignments. In `test_sample`, each row links a sample name, a FASTQ R1/R2 pair, and a library type.

Use it to confirm:

- Which FASTQ files were processed.
- Which sample each library belongs to.
- Whether each library was processed as RNA or DNA.


## `libraries/`

The `libraries/` directory contains outputs for file library pair.

Example library directories :

```text
libraries/
├── lib1_RNA_L001/
├── lib1_RNA_L002/
├── lib2_DNA_L001/
├── lib2_DNA_L002/
├── lib3_RNA_L001/
├── lib3_RNA_L002/
├── lib4_DNA_L001/
└── lib4_DNA_L002/
```

Each library directory follows this general structure:

```text
libraries/<library_id>/
├── dev/
├── metrics/
├── qc/
└── reports/
```

### `libraries/<library_id>/metrics/`

```text
metrics/
└── metrics.csv
```

These files summarize library-level metrics. In `test_sample`, these metrics focus on barcode parsing and valid full-barcode recovery.

Use these files to compare barcode performance across lanes, libraries, or samples.

### `libraries/<library_id>/qc/barcode/`

```text
qc/barcode/
├── <library_type>_barcode_occurence_heatmap.png
└── <library_type>_barcode_upset_plot.png
```

These plots show barcode recovery and barcode component completeness at the library level.

Use them to assess:

- Whether barcode usage is balanced.
- Whether complete barcode structures were recovered.
- Whether a specific lane or library has barcode dropout or barcode parsing problems.

### `libraries/<library_id>/reports/`

```text
reports/
├── <library_id>_R1_001.fastq.gz_fastqc_report.html
└── <library_id>_R2_001.fastq.gz_fastqc_report.html
```

These are FastQC reports for the raw R1 and R2 FASTQ files.

Use them to review:

- Per-base sequence quality.
- Adapter content.
- GC content.
- Read duplication.
- Other raw sequencing quality indicators.

FastQC reports are especially useful when a library has low valid barcode rate, low mapping rate, or unexpected read quality problems.

### `libraries/<library_id>/dev/`

The `dev/` directory contains intermediate files used by the pipeline. These are primarily intended for debugging, support, or advanced analysis.

For libraries:

```text
dev/<library_type>/
├── <library_id>_barcodes_mapping.parquet
├── <library_id>_corrected_barcode_reads.parquet
└── <sample_id>_reads.txt
```

These files represent read-level or barcode-level intermediate data, such as barcode correction and assignment. Customers typically do not need these files for routine analysis.

## `samples/`

The `samples/` directory contains final per-sample outputs. Each sample folder combines the relevant RNA and DNA libraries for one biological sample.

Example:

```text
samples/
├── sample1/
└── sample2/
```

Each sample directory follows this general structure:

```text
samples/<sample_id>/
├── DNA/
├── RNA/
├── counts/
├── dev/
├── metrics/
└── qc/
```

## `samples/<sample_id>/counts/`

```text
counts/
└── counts.h5mu
```

`counts.h5mu` is the main combined single-cell count object. It is the most important file for downstream analysis because it stores processed count data in a multimodal format.

Use this file for:

- Loading the sample into downstream single-cell analysis tools.
- Joint analysis of RNA gene expression and DNA amplicon-derived information.
- Integrating filtered cell-level outputs across modalities.

## `samples/<sample_id>/RNA/`

```text
RNA/
├── <sample_id>.Aligned.sortedByCoord.out.bam
├── <sample_id>.Aligned.sortedByCoord.out.bam.bai
└── <sample_id>_raw_counts.parquet
```

### RNA BAM and index

The BAM file contains RNA reads aligned to the reference genome or transcriptome. The `.bai` file is the BAM index, which allows tools to access genomic regions efficiently.

Use these files for:

- Inspecting RNA alignments in genome browsers.
- Debugging mapping behavior.
- Advanced downstream analyses that require aligned reads.

### `<sample_id>_raw_counts.parquet`

This Parquet file contains raw RNA count data before final filtering or multimodal packaging into `counts.h5mu`.

Use it when:

- A tabular count representation is needed.
- Debugging raw gene expression counts.
- Comparing raw counts to filtered counts in `counts.h5mu`.

## `samples/<sample_id>/DNA/`

```text
DNA/
├── aligned.sorted.bam
├── aligned.sorted.bam.bai
└── filtered_amplicon_counts.parquet
```

### DNA BAM and index

The BAM file contains DNA amplicon reads aligned to the reference. The `.bai` file is the corresponding index.

Use these files for:

- Inspecting amplicon read alignments.
- Reviewing target coverage.
- Debugging unexpected mapping patterns.


## `samples/<sample_id>/metrics/`

```text
metrics/
├── amplicon/
├── carryover/
├── gene_expression/
├── snp/
├── metrics.csv
└── metrics.yaml
```

### `metrics.csv``

These files summarize final per-sample metrics. The CSV is convenient for spreadsheets and scripts. The YAML is convenient for structured inspection by sample, library type, and modality.

### `metrics/amplicon/`

```text
amplicon/
├── amplicon_upset_data.csv
├── <sample_id>_anchors.csv
└── <sample_id>_mapping_rates.txt
```

These files support DNA amplicon QC.

| File | Represents |
| --- | --- |
| `amplicon_upset_data.csv` | Counts and fractions for DNA reads grouped by barcode, anchor, mapping, and high-quality mapping status. |
| `<sample_id>_anchors.csv` | Anchor-related amplicon information used to evaluate expected amplicon structure. |
| `<sample_id>_mapping_rates.txt` | Text summary of amplicon read mapping rates. |

### `metrics/gene_expression/`

```text
gene_expression/
├── <sample_id>.Log.final.out
└── <sample_id>_filtered_saturation_metrics.csv
```

These files support RNA gene expression QC.

| File | Represents |
| --- | --- |
| `<sample_id>.Log.final.out` | STAR alignment summary for the RNA alignment step. |
| `<sample_id>_filtered_saturation_metrics.csv` | Downsampled sequencing saturation and per-cell gene/UMI metrics. |

### `metrics/carryover/`

```text
carryover/
├── carryover_upset_data.csv
├── positions.csv
└── <sample_id>_mapping_rates.txt
```

These files support carryover analysis and control-region review.

| File | Represents |
| --- | --- |
| `carryover_upset_data.csv` | Counts and fractions for reads grouped by barcode, UMI, gene, mapping, cell-assignment, and carryover-related categories. |
| `positions.csv` | Positions or regions used in carryover alignment summaries. |
| `<sample_id>_mapping_rates.txt` | Text summary of mapping rates for carryover-related reads. |

### `metrics/snp/`

```text
snp/
├── <sample_id>.bcftools_stats.txt
├── <sample_id>.vcf
└── <sample_id>_snps.csv
```

These files summarize SNP or variant-related calls from the DNA analysis.

| File | Represents |
| --- | --- |
| `<sample_id>.vcf` | Variant calls in VCF format. |
| `<sample_id>_snps.csv` | Tabular summary of SNP calls or SNP counts. |
| `<sample_id>.bcftools_stats.txt` | Variant-calling statistics from bcftools. |

## `samples/<sample_id>/qc/`

```text
qc/
├── amplicon/
├── carryover/
└── gene_expression/
```

The `qc/` directory contains per-sample plots. These plots are described in more detail in the QC plot documentation, but at a high level:

| Directory | Represents |
| --- | --- |
| `qc/amplicon/` | DNA amplicon quality, mapping distribution, cell calling, and feature recovery. |
| `qc/gene_expression/` | RNA gene expression quality, cell calling, mapped read distribution, and sequencing saturation. |
| `qc/carryover/` | Carryover-related summaries and control-region signal. |

Use these plots when investigating low cell counts, low usable read fractions, weak RNA or DNA signal, unexpected mapping behavior, or possible carryover.

## `samples/<sample_id>/dev/`

```text
dev/
├── filtered_amplicon_reads.parquet
├── full_read_DNA.parquet
├── full_read_RNA.parquet
├── raw_snps.parquet
└── snps_count.parquet
```

The `dev/` directory contains intermediate or read-level files generated during processing. These are useful for support, debugging, or custom analysis, but they are not usually the first files customers need.

| File | Represents |
| --- | --- |
| `full_read_RNA.parquet` | Read-level RNA processing information. |
| `full_read_DNA.parquet` | Read-level DNA processing information. |
| `filtered_amplicon_reads.parquet` | Filtered amplicon read-level information. |
| `raw_snps.parquet` | Intermediate SNP or variant information before final summarization. |
| `snps_count.parquet` | SNP count information used to build final SNP outputs. |

## Which files should customers use?

For routine review:

- Open `experiment/experiment.html`.
- Review `experiment/metrics.csv`.
- Review `samples/<sample_id>/metrics/metrics.csv`.
- Review plots in `samples/<sample_id>/qc/`.
- Review barcode plots in `libraries/<library_id>/qc/barcode/` if barcode performance is a concern.

For downstream analysis:

- Use `samples/<sample_id>/counts/counts.h5mu` as the main multimodal count object.
- Use `samples/<sample_id>/RNA/<sample_id>_raw_counts.parquet` for raw RNA counts if needed.
- Use `samples/<sample_id>/DNA/filtered_amplicon_counts.parquet` for filtered DNA amplicon counts if needed.
- Use `samples/<sample_id>/metrics/snp/<sample_id>.vcf` or `<sample_id>_snps.csv` for variant-related review.

For troubleshooting with support:

- Share `experiment/experiment.html`.
- Share `experiment/metrics.csv`.
- Share the affected `samples/<sample_id>/metrics/` directory.
- Share the affected `samples/<sample_id>/qc/` directory.
- Share relevant `libraries/<library_id>/metrics/`, `libraries/<library_id>/qc/barcode/`, and `libraries/<library_id>/reports/` files if the issue appears library-specific.

## Notes

- The exact set of files may vary by pipeline version, assay configuration, and enabled modules.
- Files under `dev/` are primarily intermediate outputs and may change more often than customer-facing reports, metrics, counts, and QC plots.
- Empty or missing metric values can occur when a metric is not applicable or when the sample does not have enough data to estimate it reliably.
