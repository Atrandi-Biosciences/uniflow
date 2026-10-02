<p align="center">
  <img src="assets/branding/uniflow-docs-header.svg" alt="Uniflow Documentation" width="100%">
</p>

[← Repository](../README.md) · [Getting started](getting-started.md) · [Installation](installation.md) · [Inputs](inputs.md) · [Running](running.md) · [Outputs](outputs.md) · [Reports](reports.md) · [Genotyping](genotype.md) · [Support](support.md)

---

# Interpreting the QC reports

Uniflow generates a run-level `experiment/experiment.html` at run-level, and a `samples/<sample_name>/report.html` QC report for each contained sample. The table below explains the contained sections. Note that certain sections will only appear if the corresponding library, modality, or plot was produced.

## Experiment report

Open `experiment/experiment.html` first to compare libraries and samples across the run.

| Section | What it shows | How to read it |
| --- | --- | --- |
| **RNA Barcode Metrics** | Total RNA reads and the fraction with a valid full barcode. | A low valid-barcode fraction points to barcode/chemistry/read-structure loss upstream of expression quantification. |
| **RNA Gene Expression Metrics** | Raw reads, called/filtered cells, and usable-read fractions. | Use these to assess how much RNA data survives mapping, barcode/cell filtering, and UMI/gene assignment. `Filtered Cells` is the pipeline cell call from automic cell calling, unless `--force_cells` overrides this. Also consult the sample-level knee curve plot. |
| **DNA Barcode Metrics** | Total DNA reads and the fraction with a valid full barcode. | Interpret similarly to RNA barcode metrics: this primarily checks whether reads can be assigned through the expected barcode structure and white list. |
| **DNA Amplicon Metrics** | Total reads, on-target/anchored `Filtered Reads`, called cells, and `Fraction Informative Reads`. | `Filtered Reads` are reads retained after the DNA amplicon anchoring and mapping quality criteria, and feed amplicon-based cell calling. Low informative fraction can reflect poor alignment/on-target recovery. |
| **DNA Variants Metrics** | Number of discovered variant sites, agreement between callers, sites with cell support, and SNV/indel called-versus-attributed counts. | `Called` refers to the pooled discovery catalog. `Attributed` means the allele reached the per-cell layer. A difference between them means the pooled call did not obtain per-cell supporting reads under the attribution/pileup rules. |
| **Bcftools – Variant Substitution Types** | Counts of substitution classes in the called variants. | Useful as a descriptive QC view. |
| **Samtools – Flagstat** | Standard amplicon alignment/read-category counts. | Use this for a broad alignment sanity check and to spot unexpectedly unmapped, secondary, supplementary, or duplicate-heavy data where applicable. |
| **Samtools – Coverage: global stats** | Summary depth/coverage statistics across the DNA alignment. | Use for overall depth context. Aggregate depth can hide individual weak amplicons, so inspect the per-region table as well. |
| **Samtools – Coverage: stats per region** | Coverage statistics for individual reference regions/amplicons. | Look for uneven or failed amplicons that are obscured by the global average. |
| **STAR – Summary Statistics** | RNA alignment and mapping summary produced by STAR. | Use this to assess whether RNA reads map as expected before interpreting expression-derived metrics. |
| **STAR – Alignment Scores** | Distribution/summary of RNA alignment scores. | Unexpected shifts can indicate poorer read/reference agreement or alignment quality; interpret together with the STAR summary and RNA usable-read metrics. |

## Sample report

Open `samples/<sample_name>/report.html` after the experiment overview to inspect a given sample in more detail.

| Section | What it shows | How to read it |
| --- | --- | --- |
| **DNA Metrics** | Sample-level DNA Amplicon metrics exposed to the single-sample report, including total and filtered reads. | Use this as the sample-specific counterpart to the experiment table. |
| **RNA Metrics** | Sample-level RNA Gene Expression metrics exposed to the single-sample report, including total raw reads. | Use with the saturation and knee plots below. |
| **Filtered saturation** | RNA subsampling curves for median UMI counts per cell, median genes per cell, and library saturation versus mean filtered (**not** libary raw) reads per cell.  | A curve that is flattening indicates diminishing return from additional sequencing for that metric. A curve still rising suggests more depth could recover additional UMIs/genes. There is no universal pass/fail shape; interpret with the assay, cell number, and biological material. |
| **Amplicon knee plot** | Barcodes ranked by DNA amplicon read count, with amplicon-called cell status and the cell-call boundary overlaid. | A visible change in the rank curve can support the separation of high-count cells from low-count/background barcodes. A diffuse transition can occur with low-input or heterogeneous samples, so use it as a diagnostic rather than a stand-alone cell-calling rule. |
| **Gene expression knee plot** | Barcodes ranked by RNA UMI/count signal, with RNA-called cell status and the cell-call boundary overlaid. | Read it like the DNA knee plot, but for RNA expression signal. Forced cell counts, low RNA content, or substantial background can change the curve shape. For most intact cell types, a steep "cliff" is expected to separate cells from ambient background. |
| **Samtools – Coverage: global stats** | Overall DNA coverage statistics for this sample. | Gives a compact depth summary for the sample. |
| **Samtools – Coverage: stats per region** | DNA coverage by individual reference region/amplicon. | Use this to identify specific amplicons with unusually low or high coverage. |
| **Read depth for variant calling** | Distribution of variant-eligible reads per cell-amplicon pair in called cells as sequenced, on a log scale, with the variant-calling cap marked (`--downsample_reads_per_cell_amplicon`, 100 by default) and the share of pairs that reach it. | Variant calling uses at most the cap per cell and amplicon, so deeply sequenced cells do not outweigh the rest. Pairs below the cap keep all their reads. Pairs that reach the cap are used at full depth. Shares are of pairs with at least one read. Pairs with no reads are not counted. The amplicon count matrix is not capped, so its counts will not match. |

The exact sample-report sections can differ by product and by the data successfully produced for that sample. For machine-readable values, use `samples/<sample_name>/metrics/metrics.csv`; for genotype interpretation, see [Understand genotype outputs](genotype.md).

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