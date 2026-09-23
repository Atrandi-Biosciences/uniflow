<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Install Uniflow

## Requirements

Uniflow requires:

- Linux
- at least 16 CPU cores and 256 GB memory per compute node
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 25.04.8
- one supported container runtime: [Docker](https://docs.docker.com/engine/install/), [Apptainer](https://apptainer.org/docs/admin/main/installation.html), or [Singularity](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)
- internet access from the compute environment

macOS and Windows are not supported. Conda execution is not supported.

If the compute environment cannot access the internet, contact Atrandi support for an offline setup workaround.

For cluster-specific configuration examples, see [nf-core/configs](https://nf-co.re/configs/). Configuration of the STAR genome index and cluster infrastructure is the user's responsibility.

## Download a release

Downloading a release is the recommended installation method.

1. Open the [Uniflow releases page](https://github.com/Atrandi-Biosciences/uniflow/releases).
2. Download the source tarball for the latest release.
3. Extract the tarball in a location where you have write permission.
4. Open a terminal in the extracted directory.

## Clone the repository

To use the current version from the `main` branch:

```bash
git clone https://github.com/Atrandi-Biosciences/uniflow.git
cd uniflow
```

## Verify Nextflow

```bash
NXF_VER=25.04.8 nextflow -version
```

The command must report Nextflow version 25.04.8.

## Test the installation

The repository includes a sample sheet with a heavily downsampled public test dataset and public references. It verifies the installation and compute environment; it is not representative of full experimental data.

```bash
NXF_VER=25.04.8 nextflow run main.nf \
  -profile docker \
  --input_csv assets/test_samplesheet.csv \
  --product_id B4_24_v1 \
  --star_index s3://atrandi-public/public/alpha2/star_index/ \
  --amplicon_fasta s3://atrandi-public/public/alpha2/test_ref.fasta \
  --outdir <output_directory> \
  -resume
```

Replace `docker` with `apptainer` or `singularity` when required by the compute environment.


## Support

If you need help installing the pipeline, you can always reach us via our email: support@atrandi.com