<p align="center">
  <img src="docs/assets/branding/uniflow-header.svg" alt="Atrandi Biosciences — Uniflow" width="100%">
</p>

# Uniflow

Uniflow is Atrandi's Nextflow pipeline for processing data generated using Atrandi's Single-cell (RNA + DNA) Co-seq Kit. It processes paired-end FASTQ data from both RNA and DNA amplicon libraries, performs barcode correction and optional demultiplexing, runs modality-specific analysis, and produces experiment- and sample-level reports, count data, metrics, and aligned reads.

> [!WARNING]
> Uniflow is beta software for research use only. Outputs must be validated for the intended assay and sample type. Variant calling is experimental and under active development.

## Documentation

- [Getting started](docs/getting-started.md)
- [Installation](docs/installation.md)
- [Prepare inputs](docs/inputs.md)
- [Run Uniflow](docs/running.md)
- [Understand outputs](docs/outputs.md)
- [Output definitions](docs/specs/outputs.md)
- [QC report files](docs/reports.md)
- [Genotyping](docs/genotype.md)
- [Troubleshooting and support](docs/support.md)

## Issues and support

Report reproducible code issues through [GitHub Issues](https://github.com/Atrandi-Biosciences/uniflow/issues).

Uniflow is distributed under the terms in [LICENSE.md](LICENSE.md).
