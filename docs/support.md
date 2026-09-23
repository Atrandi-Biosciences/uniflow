<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Troubleshooting and support

Pipeline users should contact their Atrandi Field Application Scientist (FAS). The files and information to provide depend on whether the problem concerns assay performance or pipeline execution.

## Sample or assay issue

Use this route when the pipeline completes but a sample or modality fails, produces unexpected results, or shows poor assay performance.

Send your FAS contact:

- `experiment/qc_report.tar.gz`
- the affected sample or library identifier
- a short description of the problem

The support archive contains assay-performance metrics, reports, run parameters, and QC plots. It does not contain reads or information intended to reveal the underlying biology.

The archive does contain sample names and library IDs. Use neutral identifiers in the input sample sheet if these names are sensitive.

## Pipeline issue

Use this route when the pipeline stops before completing. Send the `.nextflow.log` file from the directory where Uniflow was launched to your FAS contact.

Files in `pipeline_info/` contain additional execution diagnostics. Provide them if your FAS contact requests them; they are not required for the initial report.

## Developer and code issues

Developers can report reproducible code defects or request changes through [GitHub Issues](https://github.com/Atrandi-Biosciences/uniflow/issues).

Do not use public GitHub issues for customer support or attach files containing sensitive identifiers.

## General support channel

You can always contact us via support@atrandi.com for support regarding your assay or pipeline issues.