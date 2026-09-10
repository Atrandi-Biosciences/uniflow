# Chemistry Specification

## Overview

A chemistry definition describes the structure of sequencing reads for a specific assay. It provides a declarative configuration that allows the pipeline to parse reads, identify functional sequence elements, and associate them with external reference libraries without requiring chemistry-specific code.

The design separates the physical layout of a read from its biological interpretation, making new chemistries easy to add and maintain.

---

## Design Principles

The chemistry specification is designed to be:

- **Declarative**: describes read structure without encoding parsing logic.
- **Readable**: symbolic patterns provide a compact representation.
- **Modular**: layout, biological meaning, references, and operational behavior are separated.
- **Extensible**: new chemistries can be introduced through configuration alone.
- **Versioned**: multiple assay designs can coexist within the same pipeline.

---

# Chemistry Definition

Each chemistry is identified by a unique, versioned name.

```yaml
B4-24_v1:
```

A typical naming convention is:

```
<Design-shortname>_v<Version>
```

This allows multiple chemistry versions to coexist while preserving backward compatibility.

An optional `description` field provides a human-readable summary.

```yaml
description: "24 barcodes, 4 segments, version 1"
```

---

## Read Structure

The `structure` section defines one or more read library types.

```yaml
structure:
  rna:
  dna:
```

Each read type can have an independent organization while sharing common concepts such as barcodes and linkers.

### Offset

```yaml
offset: 0
```

The offset specifies the starting position for parsing. Bases before the offset are ignored.

### Pattern

The `pattern` defines the read layout using symbolic characters.

Example:

```
DDDDDDDDOOOOCCCCCCCCNNNNBBBBBBBBMMMMAAAAAAAALLLLHHHHHHHHHHHHHHUUUUUUUU
```

Each character represents a nucleotide position, while consecutive identical characters form a sequence segment.

The parser determines segment lengths directly from the pattern, making the specification compact and easy to modify.

---

## Feature Mappings

Pattern symbols are assigned biological meaning through the `mappings` section.

Example:

```yaml
mappings:
  barcode:
    D: Barcode D
    C: Barcode C
    B: Barcode B
    A: Barcode A
```

Mappings group sequence segments into logical feature classes, such as:

- **Barcode**: sample or molecule identifiers.
- **Linker**: spacer sequences between barcodes.
- **Handle**: capture or amplification sequences.
- **UMI**: unique molecular identifiers.

Separating the pattern from its biological interpretation allows the parser to remain generic and extensible.

`barcode` and `handle` are special features which will be corrected if a reference is provided.
`UMI` is corrected by default.

---

## Reference Libraries

Feature classes can use external sequence references.

```yaml
references:
  barcode:
    path: ...
  handle:
    path: ...
```

Reference files contain valid sequence definitions used for matching and error correction. Those files can be `parquet`, `csv`, `csv.gz`. We require parquet files for larger references (>500KB compressed csv)

### Reference Headers

`barcode`:
- `symbol`: Letter used in the pattern string
- `sequence`: Reference sequence
- `plate_coordinate` (optional): Coordinate of the well. This is only used for barcode refs in plates. Example: A1

`handle`:
- `symbol`: Letter used in the pattern string
- `sequence`: Reference sequence

---

## Demultiplexing

The `demultiplexing` section specifies which feature components is used for sample assignment.

```yaml
demultiplexing:
  parts: A
```

The `parts` field refers to which symbol the sample demultiplexing will occur.


