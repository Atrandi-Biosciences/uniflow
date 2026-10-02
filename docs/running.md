<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Reports](reports.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Run Uniflow

Run Uniflow from the extracted release or cloned repository directory:

```bash
NXF_VER=25.04.8 nextflow run main.nf \
  -profile docker \
  --input_csv <sample_sheet.csv> \
  --product_id <product_id> \
  --star_index <star_index> \
  --amplicon_fasta <amplicon_reference.fasta> \
  --outdir <output_directory> \
  -resume
```

Replace `docker` with `apptainer` or `singularity` to use another supported container runtime.

The input and output paths may be local paths or S3 URLs. See [Prepare inputs](inputs.md) for product IDs, sample-sheet requirements, and reference behavior.

Runtime depends on sample size, sequencing depth, and the available compute resources and configuration. With the minimum recommended system resources, a typical RNA + DNA Co-Sequencing sample is expected to complete in approximately 8–10 hours. This estimate assumes approximately 20,000 cells, 20,000 RNA reads per cell, and up to 20 amplicons, each sequenced to up to 1,500 reads per cell.

## Run a single library type

Most experiments provide both RNA and DNA libraries with both references. Uniflow can still process a sample when one library type failed: remove the failed library rows from the sample sheet. If no libraries of that type remain anywhere in the sheet, its corresponding reference argument may also be omitted. Uniflow will process the remaining library type normally.

Always use `-resume` when restarting a run. Nextflow will reuse completed work when the inputs and relevant parameters have not changed.

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
