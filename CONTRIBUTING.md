# Contributing to uniflow

This document is underlying naming conventions, code structure and pipeline architecture to simplify contributing to any developer.

Hopefully this document will help you know exactly where a new piece of code should be and help you add a new feature in no time.

# Key concepts

The uniflow pipelines stands on a few core concepts that help define its structure.

1. A library. A library is a collection of reads coming from the same physical sequencing library (same index).
2. A sample. A collection of reads that originate from the same conditions (biological, experiment and/or technical)
3. A library type. A type of library, arbitrary name to distinguish different molecular origins
4. A modality. A type of obervation originating from a process applied to a library. A common well known one is `gene expression`


# Subworkflows

Subworkflows are the first entry point coming from main.
The goal is to keep those entrypoints independent of other potential entrypoints keeping the tree of dependencies clean.
Example:

`dna.nf` contains all code paths related to DNA

## Modality entrypoints

Keep modality-level orchestration in `subworkflows/<modality>.nf`.
These files should read like the table of contents for a modality: include the
high-level subworkflows, pass shared inputs to them, connect only the outputs
that are needed across those subworkflows, and normalize the final emits.

For example, `subworkflows/rna.nf` owns RNA-level routing. It includes
`gene_expression` and `cross_processing`, takes the RNA libraries plus shared
resources such as the STAR index, whitelist, sample barcode mappings, trimming
length, amplicon fasta, and run switches, then emits grouped RNA outputs.

When adding a new modality or a new high-level input:

1. Add product-level modality switches and global resources in `main.nf`.
2. Pass only high-level inputs into the modality entrypoint.
3. Put modality-specific feature logic in an included subworkflow directory,
   such as `subworkflows/gene_expression/`.
4. Keep process-level implementation in `modules/`.
5. Emit stable grouped channels from the modality entrypoint.

Use this shape for modality emits:

1. `metrics`, `qc`, and `dev`: `[meta, modality, [files]]`, grouped by
   `[meta, modality]`.
2. `out`: `[meta, [files]]`, grouped by `meta`.
3. `metrics_parquet`: `[meta, [files]]`, grouped by `meta`.
4. Modality-specific primary outputs may keep explicit names, for example
   `counts_h5ad`.

If a feature subworkflow produces files for another feature or modality, name
those emits for the consumer first, for example
`cross_processing_merged_reads`. This keeps cross-feature coupling visible at
the boundary instead of hidden inside process code.

# Modules

# Processes

## Naming Conventations

We typicall recommend using a small descriptive process name such as starting by a verb at the present tense followed by the subject. Example: `CORRECT_BARCODE`

# Metrics handling and propagation

# File formats

Here we will recommend file formats for commonly used temporary files or outputs. You are free to introduce new ones but it must be accepted by Atrandi ahead of time and included here.

1. parquet: Bread and buther file for storing temporary files. They are compressed tabular files that allow to use polars or duckdb.
2. csv: For small files or files that need to be directly accessed by a person without any complicated interface. Preferred for metrics outputs
3. bam: Generally accepted file format for aligned reads
4. vcf: WIP
5. html: Used for interactive reports
6. png: Used for any image that needs to be outputed by the pipeline.

# [WIP] Output files and formats

Here is the general structure of a pipeline run

```
├── experiment: All files reporting an aggregation event of many samples
├── libraries: Outputs related to libraries, each subfolder is a library
├── metadata: Files related to the run, params and samplesheet
└── samples: Outputs related to samples, each subfolder is a sample
```



# Container image (Wave)

The Python runtime shared by most processes is a single container image built
with Seqera Wave. Its full, version-locked environment is in
[`bin/conda.yaml`](bin/conda.yaml): the conda-forge/bioconda channels, a pinned
`python`/`pip`, and an explicit pin for every pip dependency. The same file is
both the Wave build recipe and the local dev environment spec, so a contributor
can reproduce the image contents with
`micromamba env create -f bin/conda.yaml` or an equivalent command.

Each process references the image by its Wave name, for example:

```
container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
```

The tag is a content hash assigned by Wave, so it changes whenever the
environment changes; every process module on this image must point at the same
tag. The UpSet-plot QC processes are the exception, running on a separate
`uniflow-qc` image (`bin/conda_qc.yaml`) that adds upsetplot pinned to pandas <3.

## Building and bumping the image

We build manually through the Wave website for now (no Wave CLI yet). To change
a dependency and use a new image, follow these three steps:

1. Bump the version. Edit `bin/conda.yaml`: add the new package or change an
   existing pin, keeping every dependency explicitly version-locked (no unpinned
   entries). Recreate your local env to confirm it resolves:
   `micromamba env create -f bin/conda.yaml`.
2. Post the env file on Wave. Open the builder at https://wave.seqera.io, choose
   the edit environment option, and paste the full contents of `bin/conda.yaml`.
   Build. Wave returns an image reference of the form
   `community.wave.seqera.io/library/uniflow:<new-tag>`.
3. Update the images. Replace the old tag with `<new-tag>` in every process
   `container` directive under `modules/`; they must all match. Find them with:

   ```bash
   grep -rn "community.wave.seqera.io/library/uniflow:" modules
   ```

   or rewrite them in one pass (macOS `sed` for example). If `xargs` is not available,
   you can install `build-essential` and `findutils` in your favorite *nix OS.

   ```bash
   # example for macOS
   grep -rl "community.wave.seqera.io/library/uniflow:" modules \
     | xargs sed -i '' 's|uniflow:<old-tag>|uniflow:<new-tag>|g'
   ```

On pull requests, CI runs a "Container image manifest check" that resolves every
image referenced by the pipeline (via `skopeo`), so a mistyped or unbuilt tag
fails the build before merge.

## Conventional Commits
Conventional Commits is a lightweight standard for writing human- and machine-readable Git commit messages.
- https://www.conventionalcommits.org/en/v1.0.0/
- https://python-semantic-release.readthedocs.io/en/latest/ is used for the technical implementation.

Python-semantic-release process the commits and calculates the new version.
This version is then updated in the nextflow.config file. New features and bug fixes are written to the CHANGELOG.md file.
All modified files are pushed to the main branch, and a tag with the new version is set.

* * * Important!
Python-semantic-release cannot read the title of a pull request. When merging with the squash method, it is important to create the merge commit according to the conventional commit rules. Examples are described below.

## PR-Titel and Version (after Merge to main)

- `feat!:` / `fix!:` → Major (x.0.0)
- `feat:`  → Minor (1.x.0)
- `fix:`   → Patch (1.0.x)

## Allowed Types

| Type     | Description                      | Version Impact (Later) |
|----------|----------------------------------|------------------------|
| feat     | New feature                      | Minor                  |
| feat!    | New feature with breaking change | Major                  |
|          | (not backwards compatible)       |                        |
| fix      | Bug fix                          | Patch                  |
| fix!     | Bug fix with breaking change     | Major                  |
|          | (not backwards compatible)       |                        |
| docs     | Documentation                    | –                      |
| test     | Tests (e.g., nf-test) / CI-CD    | –                      |
| refactor | Code change with no behavior fix | –                      |
| revert   | Revert                           | –                      |

## Examples
   type'optional'(...): commit-message
   feat: add end-to-end nf-test for mixed sample sheet
   feat(parser): add end-to-end nf-test for mixed sample sheet
   feat!: rename parameter --reads to --input
   fix: correct barcode_pool_size handling in demux
   test: add snapshot for variant calling output / add commitlint step to pull requests
   docs: document required Bitbucket repository variables
   refactor: extract qc processes into subworkflow
   
## BITBUCKET_TOKEN
python-semantic-release needs the token to commit the new version to the main branch and set a tag.
- Type: Repository Access Token
- Created at: Repository settings → Access tokens
- Stored at: Repository settings → Pipelines → Repository variables (secured)
- Scopes: Pull requests Read
- optional: Repositories Write (for Git tags / version commits)
- Purpose: Read the PR title via API after merging into main to determine major/minor/patch