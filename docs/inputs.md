<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Prepare inputs

Uniflow supports RNA-only, DNA-only, and combined RNA + DNA Co-seq experiments.

Input and output paths may use either the local filesystem or S3 URLs.


## References

| Libraries | Required reference |
| --- | --- |
| RNA | STAR genome index (`--star_index`) |
| DNA | Amplicon FASTA (`--amplicon_fasta`) |

Uniflow processes a library type only when both its FASTQ inputs and required reference are provided. If either is missing, that library type is skipped. For example, an experiment with RNA and DNA libraries but only a STAR index runs the RNA analysis and skips DNA analysis.

Customers are responsible for creating and validating their STAR index. The index must be built with STAR 2.7.11b, matching the version used by Uniflow.

## Create the amplicon reference

Provide one FASTA record per amplicon. Each sequence must contain the complete expected amplicon inlcuding the forward primer, target region, and reverse primer sequence.

```text
>AMPLICON_1
ACGTACGT[...]ACGTACGT
```

Amplicon IDs must be unique and may contain letters, numbers, and underscores. Sequences may contain `A`, `C`, `G` or `T`; uppercase is recommended.

## Create the sample sheet

Provide a CSV file to `--input_csv`. Each row represents one paired-end FASTQ input. Reuse a `sample_name` across rows to aggregate multiple libraries or sequencing lanes into one sample.

Start from the [test sample sheet](../assets/test_samplesheet.csv) or create a CSV with these required columns:

| Column | Description |
| --- | --- |
| `sample_name` | Sample to which the FASTQ pair belongs. May contain letters, numbers, underscores, and hyphens. |
| `fastq_1` | Path to the R1 FASTQ file. |
| `fastq_2` | Path to the matching R2 FASTQ file. |
| `library_type` | `RNA` or `DNA`. |

A sample having **both** RNA and DNA data must hold the **same** sample_name to allow the pipeline to merge the modalities downstream.

Example:

| `sample_name` | `fastq_1` | `fastq_2` | `library_type` |
| --- | --- | --- | --- |
| `sample1` | `/data/Lib1_RNA_L001_R1_001.fastq.gz` | `/data/Lib1_RNA_L001_R2_001.fastq.gz` | `RNA` |
| `sample1` | `/data/Lib2_DNA_L001_R1_001.fastq.gz` | `/data/Lib2_DNA_L001_R2_001.fastq.gz` | `DNA` |

FASTQ files must be gzip-compressed and named in the following format (Illumina standard naming convention):

```text
<library_id>_L00X_R1_001.fastq.gz
<library_id>_L00X_R2_001.fastq.gz
```

The R1 and R2 files in each row must have the same library ID and lane.

### Optional columns

| Column | Description |
| --- | --- |
| `demultiplexing_indices` | Semicolon-separated plate coordinates, such as `A1;A2;A3`, used to split a library into samples. Omit the column when demultiplexing is not needed. RNA and DNA rows belonging to the same sample must use exactly the same coordinates. |
| `force_cells` | Positive integer overriding automatic cell thresholding. Use only when the automatic threshold is unsatisfactory. |


## Select a product

Pass one of these values to `--product_id`:

| Product ID | Plate barcodes | Variant calling |
| --- | ---: | --- |
| `B4_24_v1` | 24 | Enabled |
| `B4_96_v1` | 96 | Enabled |
| `B4_24_v1-no-variant` | 24 | Disabled |
| `B4_96_v1-no-variant` | 96 | Disabled |

Use a `-no-variant` product when the data is not suitable for variant calling. The remaining `gene expression` and `amplicon` analysis still runs.
