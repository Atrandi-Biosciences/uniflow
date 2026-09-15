# Uniflow

CONFIDENTIAL – EARLY ACCESS SOFTWARE 

NO REDISTRIBUTION 

Do not copy, publish, upload, sublicense, transfer, or distribute this software to 

third parties without prior consent by Atrandi Biosciences. 
Uniflow is Atrandi’s RNA+DNA co-sequencing pipeline built with Nextflow. It processes paired-end FASTQs from RNA (cDNA) and DNA (amplicon) libraries, performs barcode correction and demultiplexing, runs modality-specific alignment/quantification, and produces per-sample outputs, QC, and aggregated metrics. 

# DISCLAIMER

This is an alpha version of the pipeline to allow early adopters to gain insight into their data. This product is for research use only. There might be bugs, this software is provided as is with no guarantee of results. 

## Variant calling

Uniflow runs **germline SNV** calling on amplicon BAMs via multiple callers in parallel (freebayes, bcftools, VarDict, GATK HaplotypeCaller by default; LoFreq is wired and opt-in via `--variant_callers`). Per-caller VCFs are normalized to SNVs only, merged into a single sample-level consensus VCF carrying `INFO/NCALLERS` and per-caller `FORMAT/{VAF,AD,DP,GT}`, and attributed per cell via pysam pileups. The MultiQC report surfaces per-caller and consensus metrics in two tables plus a caller-agreement UpSet plot per sample.

Somatic mode is out of scope for v1.

**Breaking change vs. earlier alphas:** the legacy pseudo-bulk `bcftools` SNP path (`PILEUP_READS` → `snps.csv`) was removed. Pipelines that consumed `snp_calling/*.h5ad` should switch to the per-caller outputs under `samples/{sample}/DNA/variants/{caller}/per_cell/`.

The current SNP caller implementation is currently only supporting SNP variants. We use a pseudo bulk approach on the whole sample to call potential variants. We then use those to find the single cell variants. We are actively working on a single-cell caller for our beta phase.

# Requirements

## Environments

The pipeline is compatible with both [docker](https://docs.docker.com/engine/install/) and [apptainer](https://apptainer.org/). Before trying to run the pipeline, make sure the person responsible for data processing is aware of what backend is available.

**The conda options is not enabled for this pipeline**, docker or apptainer is **required**.

## Platform

This pipeline has only been tested on Linux. MacOS, Windows or other operating systems are **not** supported by Atrandi. 

## System

The current minimum requirement for a node is `16 cores` and `256GB` of memory.

If you need help with setting up your cluster environment, you can follow the [official nextflow guidelines](https://www.nextflow.io/docs/latest/executor.html)

[Here](https://nf-co.re/configs/) are a list of predefined configuration from different institutions.

see also [configuration](#configuration)

## Internet connection

Uniflow currently requires an internet connection to download containers to run.

The current SNP caller implementation is quite poor. We use a pseudo bulk approach on the whole sample to call potential variants. We then use those to find the single cell variants. We are actively working on a single-cell caller for our beta phase.

# Requirements

## Installation
To use our pipeline you will need to first:

1. [install nextflow](https://nextflow.io/docs/latest/install.html#install-page).
2. Install [docker](https://docs.docker.com/engine/install/) or [apptainer](https://apptainer.org/) if not already done.
3. Download a fasta and GTF file for your species.
4. Generate a reference with [STAR](https://github.com/alexdobin/STAR) (v2.7.11b)
5. Manually generate a fasta file for your amplicon reference using this [guideline](#amplicon-reference)

**The conda options is not enabled for this pipeline**, docker or apptainer is **required**.

## Configuration

There is a configuration file at the root of the pipeline called `nextflow.config`. This is where you might do modifications to get the pipeline compatible with your local environment.

You need to set the `MAX_CPU` available on your compute in the `resources.config` file.


## Samplesheet

The pipeline requires you to provide a samplesheet with the following columns:

- `sample_name`: Defines the sample name. A sample could be any condition (biological or technical). Samples names will be included in metrics outputs that might be requested by our support team, please don't use descriptive sample names.
  - validation: Sample names need to start by a letter and they **cannot** contain spaces. Please make use of `_`. Example: `sample_1`
  - **IMPORTANT**: Do not use confidential naming patterns as these names will be in QC plots potentially share with customer support.
- `library_type`: `RNA` or `DNA`. RNA for cDNA and DNA for amplicon. This will define what your library files originated from.
  - validation: File names including `R1`, `R2`, or `L00X` strings will be rejected by the pipeline.
- `fastq_1`: Path to the R1 fastq.gz file of a specific illumina index. The R1 file contains sequences to align to a reference.
  - **IMPORTANT**: Do not use confidential naming patterns as these names will be in QC plots potentially share with customer support.
- `fastq_2`: Path to the R2 fastq.gz file of a specific illumina index. The R2 file contains sequences for the SPC barcodes.
  - **IMPORTANT**: Do not use confidential naming patterns as these names will be in QC plots potentially share with customer support.
- `force_cells` <Optional>: Number of cells to threshold your samples on.

There is test samplesheet provided that you can use to make sure your infrastructure is ready here: `assets/test_samplesheet.csv`

**This test dataset is heavily downsampled and does not reflect real data.**

```
`NXF_VER=25.04.8 nextflow run main.nf --input_csv assets/test_samplesheet.csv --star_index s3.... --amplicon_fasta s3.... --outdir PATH_TO_OUTPUT -resume`
```


## Amplicon reference

Please make sure the fasta reference is written as follows:

```
>AMPLICON_ID
FORWARDNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNREVERSE
```

Another way to represent it is:

```
FORWARD_PRIMER_SEQUENCE + TARGET_SEQUENCE + REVERSE_PRIMER_SEQUENCE 
```

Make sure the `amplicon_id` is in all caps.
- FORWARD is the forward primer sequence
- N: sequence between the primers
- REVERSE is the reverse primer sequence

# Running the pipeline

This pipeline has only been tested on Linux. MacOs and Windows are **not** supported by Atrandi.

1. Unzip the provided pipeline at a location you have write permissions.
2. `cd uniflow`
3. Run the pipeline using the following command: `NXF_VER=25.04.8 nextflow run main.nf --input_csv SAMPLE_SHEET_PATH --star_index PATH_TO_STAR_INDEX --amplicon_fasta PATH_TO_AMPLICON_FASTA --outdir PATH_TO_OUTPUT -resume`

## Uniflow parameters

- --input_csv: Path to the csv file you prepared
- --star_index: Path to the folder that holds the STAR index files
- --amplicon_fasta: Path to the fasta amplicon reference
- --outdir: Path to the output folder where results are going to be written to

Note: All these paths can be s3 path, uniflow will download them for you.

# Outputs

Output documentation and troubleshooting guidelines can be found [here](docs/outputs.md)

# Quality control and communication with customer suppoer

You can find a detailed explanation and interpretation of all our QC plots [here](docs/plots.md)


