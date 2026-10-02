<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Reports](reports.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Understand Uniflow outputs

Start with the experiment report for a run-wide overview, then open a sample report for details about a specific sample.

For file format and field level details, see the [output definitions specification](specs/outputs.md).

The exact files returned depend on the libraries, references, and product selected. A combined RNA+DNA run with variant calling produces this customer-visible structure:

```text
<output_directory>/
├── experiment/
│   ├── experiment.html
│   ├── metrics.csv
│   └── qc_report.tar.gz
├── libraries/
│   └── <library_id>/
│       ├── metrics/metrics.csv
│       └── reports/*_fastqc_report.html
├── metadata/
│   ├── chemistry.json
│   ├── params.json
│   └── sample_sheet.csv
├── pipeline_info/
│   ├── report_<timestamp>.html
│   ├── timeline_<timestamp>.html
│   └── trace_<timestamp>.txt
└── samples/
    └── <sample_name>/
        ├── report.html
        ├── status.csv
        ├── counts/
        │   ├── counts.h5mu
        │   └── counts.zarr/
        ├── DNA/
        │   ├── aligned.sorted.bam
        │   ├── aligned.sorted.bam.bai
        │   ├── <sample_name>.consensus.vcf.gz
        │   └── <sample_name>.consensus.vcf.gz.tbi
        ├── RNA/
        │   ├── <sample_name>.Aligned.sortedByCoord.out.bam
        │   └── <sample_name>.Aligned.sortedByCoord.out.bam.bai
        ├── dev/
        │   ├── amplicon/
        │   ├── gene_expression/
        │   └── variant/
        │       ├── variants_indel.parquet
        │       └── variants_snv.parquet
        └── metrics/
```

Standalone sample and library QC directories are not published for customer runs. Relevant plots are included in the HTML reports and QC support bundle.

## Experiment overview

| Output | Purpose |
| --- | --- |
| `experiment/experiment.html` | Primary run-level report. Open this first to review all libraries and samples. |
| `experiment/metrics.csv` | Combined run-level metrics in a machine-readable table. |
| `experiment/qc_report.tar.gz` | Diagnostic bundle for Atrandi support. It contains the experiment report, metrics, run parameters, and QC plots. |

## Library outputs

Each input FASTQ pair has a library directory.

| Output | Purpose |
| --- | --- |
| `libraries/<library_id>/metrics/metrics.csv` | Barcode-processing metrics for the library. |
| `libraries/<library_id>/reports/*_fastqc_report.html` | Read-quality reports for the R1 and R2 FASTQs. |

## Sample outputs

| Output | Purpose |
| --- | --- |
| `samples/<sample_name>/report.html` | Detailed report for one sample. Open this after reviewing the experiment report. |
| `samples/<sample_name>/status.csv` | Processing status for the sample's modalities. |
| `samples/<sample_name>/metrics/metrics.csv` | Combined metrics for the sample. |
| `samples/<sample_name>/metrics/<modality>/` | Supporting metrics from the RNA or DNA analysis. |
| `samples/<sample_name>/RNA/` | Coordinate-sorted RNA BAM and BAI files when RNA was processed. |
| `samples/<sample_name>/DNA/` | DNA BAM and BAI files. When enabled, this directory also contains the consensus VCF and its index. |
| `samples/<sample_name>/dev/<modality>/` | Exploratory outputs and formats for power users. These may change between releases and are not the recommended starting point. |

## Count data

Both count formats contain the combined per-sample modalities:

| Output | Use |
| --- | --- |
| `samples/<sample_name>/counts/counts.h5mu` | Conventional MuData file for local downstream analysis. |
| `samples/<sample_name>/counts/counts.zarr/` | Zarr representation suitable for direct access from S3. |

Use the format supported by the downstream analysis environment; neither format contains additional analysis results unavailable in the other.

Each modality is stored under its own key: `gene_expression` and `amplicon`, plus `snv` and `indel` when variant calling is enabled.

## Variant and genotype data

When variant calling is enabled, each sample with DNA data also has:

| Output | Purpose |
| --- | --- |
| `samples/<sample_name>/counts/counts.h5mu`, modalities `snv` and `indel` | Per-cell genotypes, genotype qualities and read counts for SNVs and indels. Start here for genotype analysis. |
| `samples/<sample_name>/DNA/<sample_name>.consensus.vcf.gz` | Variants discovered from the pooled reads of all cells, with one sample column per variant caller. Contains no per-cell genotypes. |
| `samples/<sample_name>/dev/variant/variants_snv.parquet` | Per-cell SNV table with read counts, qualities, genotype likelihoods and caller support. |
| `samples/<sample_name>/dev/variant/variants_indel.parquet` | The same table for insertions and deletions. |
| `samples/<sample_name>/metrics/metrics.csv` | Variant discovery and attribution metrics, in the rows whose `modality` is `variants`. |

See [Understand genotype outputs](genotype.md) for how genotypes are called, what each field means, and how to avoid common misreadings.

## Metadata and execution details

`metadata/` records the resolved sample sheet, product chemistry, and run parameters needed to understand or reproduce the run.

`pipeline_info/` contains Nextflow execution diagnostics. These files describe task resource use, timing, and failures; they are primarily useful when troubleshooting a pipeline-level problem.

## Working with Uniflow outputs

For examples of loading Uniflow count data and performing downstream analysis, see the relevant Uniflow analysis [vignettes](https://github.com/Atrandi-Biosciences/vignettes).


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