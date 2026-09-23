<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Genotyping](genotype.md) · [Support](support.md)

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

## Run a single library type

Most experiments provide both RNA and DNA libraries with both references. Uniflow can still process a sample when one library type failed: remove the failed library rows from the sample sheet. If no libraries of that type remain anywhere in the sheet, its corresponding reference argument may also be omitted. Uniflow will process the remaining library type normally.

Always use `-resume` when restarting a run. Nextflow will reuse completed work when the inputs and relevant parameters have not changed.
