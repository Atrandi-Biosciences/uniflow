from collections.abc import Mapping
from enum import StrEnum
from typing import Any, Final

# Barcode naming
BARCODE_POSITION: Final[str] = "barcode_position"


class BarcodePosition(StrEnum):
    """Barcode position in sequence."""

    ONE = "1"
    TWO = "2"
    THREE = "3"
    FOUR = "4"


# Barcode name
BARCODE_D: Final[str] = "D"
BARCODE_C: Final[str] = "C"
BARCODE_B: Final[str] = "B"
BARCODE_A: Final[str] = "A"
BARCODE: Final[str] = "barcode"
LINKER_1: Final[str] = "linker_1"
LINKER_2: Final[str] = "linker_2"
LINKER_3: Final[str] = "linker_3"
LINKER_4: Final[str] = "linker_4"

# Library naming
SEQUENCE: Final[str] = "sequence"
READ_NAME: Final[str] = "read_name"

MODALITY: Final[str] = "modality"


class LibraryType(StrEnum):
    """Library types supported by the pipeline."""

    DNA = "DNA"
    RNA = "RNA"


class Modality(StrEnum):
    """Pipeline modularity levels."""

    GENE_EXPRESSION = "gene_expression"
    AMPLICON = "amplicon"
    BARCODE = "barcode"
    SNV = "snv"
    INDEL = "indel"
    VARIANTS = "variants"
    CARRYOVER = "carryover"
    NONE = "none"


class Features(StrEnum):
    """Standard feature column names."""

    GENE_ID = "gene_id"
    GENE_NAME = "gene_name"
    AMPLICON_ID = "amplicon_id"
    SNP_ID = "snp_id"


# RNA constants
UMI: Final[str] = "umi"

PRIMER: Final[str] = "primer"

# DNA constants
CHROM: Final[str] = "chrom"

COUNTS: Final[str] = "counts"
READS: Final[str] = "reads"
INDEX: Final[str] = "index"
FRACTION: Final[str] = "fraction"

MAPPING_QUALITY: Final[str] = "mapping_quality"
MAPPED: Final[str] = "mapped"
HIGH_QUALITY_MAPPED: Final[str] = "high_quality_mapped"
UNMAPPED: Final[str] = "unmapped"
HAS_UMI: Final[str] = "has_umi"
HAS_BARCODE: Final[str] = "has_barcode"
HAS_GENE: Final[str] = "has_gene"
IN_CELL: Final[str] = "in_cell"

# Sample sheet naming
SOURCE_ID: Final[str] = "source_id"
DEMULTIPLEXING_INDICES: Final[str] = "demultiplexing_indices"
FASTQ_1: Final[str] = "fastq_1"
FASTQ_2: Final[str] = "fastq_2"
LIBRARY_TYPE: Final[str] = "library_type"
SAMPLE_NAME: Final[str] = "sample_name"
LIBRARY_ID: Final[str] = "library_id"
EXPECTED_CELLS: Final[str] = "expected_cells"
PROCESS_NAME: Final[str] = "process_name"

# Barcode plotting constants
TOTAL_COUNTS: Final[str] = "total_counts"
N_FEATURES: Final[str] = "n_features"
CONFIDENCE_SCORE: Final[str] = "confidence_score"
IS_CELL: Final[str] = "is_cell"
RANK: Final[str] = "rank"

# Human readable labels registry
# TODO: Replace this by a configuration file
HUMAN_READABLE_NAMES: dict[Any, str] = {
    BARCODE_A: "Barcode A",
    BARCODE_B: "Barcode B",
    BARCODE_C: "Barcode C",
    BARCODE_D: "Barcode D",
    PRIMER: "RNA Handle",
    LibraryType.DNA: "DNA",
    LibraryType.RNA: "RNA",
    Modality.GENE_EXPRESSION: "Gene Expression",
    Modality.AMPLICON: "Amplicon",
    Modality.BARCODE: "Barcode",
    Modality.SNV: "SNV",
    Modality.INDEL: "Indel",
    Modality.VARIANTS: "Variants",
    Modality.CARRYOVER: "Carryover",
}
HUMAN_READABLE_NAMES_FROZEN: Mapping[Any, str] = dict(
    HUMAN_READABLE_NAMES
)  # shallow-immutable view


def label_for(key: Any, default: str | None = None) -> str:
    """
    Return the human-readable label for a given registry key.

    Args:
        key: Column key or enum (e.g., BARCODE_D or LibraryType.DNA).
        default: Fallback string if key is not registered. If None, returns str(key).

    Returns:
        The human-readable label for the key.
    """
    if key in HUMAN_READABLE_NAMES_FROZEN:
        return HUMAN_READABLE_NAMES_FROZEN[key]
    return default if default is not None else str(key)


# Genotype layer names
LAYER_GT: Final[str] = "GT"
LAYER_GQ: Final[str] = "GQ"
