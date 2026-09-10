"""Compile chemistry YAML definitions into runtime configuration models."""

from collections import Counter
from pathlib import Path
from typing import Any

import polars as pl

from lib.chemistry.components import (
    DemultiplexingStrategy,
    LibraryInfo,
    ReferenceInfo,
    SequencePart,
)

from lib.chemistry.models import (
    LIST_OF_FEATURE_CLASSES,
    REQUIRED_COLUMNS_BY_FEATURE_TYPE,
)
from lib.chemistry.chemistry import ChemistryDef
from lib.chemistry.common_const import (
    DESCRIPTION,
    STRUCTURE,
    OFFSET,
    PATTERN,
    MAPPINGS,
    REFERENCES,
    SYMBOL,
    DEMULTIPLEXING,
    REFERENCE_READERS,
    PARTS,
)


def compile_chemistry_definitions(
    chemistry_definitions: dict[str, Any], project_root: Path
) -> dict[str, dict[str, Any]]:
    """Compile all chemistry definitions from source YAML data.

    Args:
        chemistry_definitions: Parsed top-level chemistry YAML.
        project_root: Project root used to resolve reference paths.

    Returns:
        JSON-serializable runtime chemistry definitions.
    """
    compiled: dict[str, dict[str, Any]] = {}
    for chemistry_id, config in chemistry_definitions.items():
        chemistry_def = compile_chemistry_definition(
            chemistry_id=chemistry_id,
            config=config,
            project_root=project_root,
        )
        compiled[chemistry_id] = chemistry_def.to_dict()
    return compiled


def deduct_demultiplexing_type(
    chemistry_def: ChemistryDef,
    demultiplexing_in_libraries: int,
    demultiplexing_parts: list[str],
) -> DemultiplexingStrategy | None:
    n_libraries = len(chemistry_def.libraries)
    if n_libraries > demultiplexing_in_libraries:
        return DemultiplexingStrategy(type="partial", symbols=demultiplexing_parts)
    elif n_libraries == demultiplexing_in_libraries:
        return DemultiplexingStrategy(type="full", symbols=demultiplexing_parts)
    return None


def compile_chemistry_definition(
    chemistry_id: str, config: dict[str, Any], project_root: Path
) -> ChemistryDef:
    """Compile one source chemistry definition.

    Args:
        chemistry_id: Top-level chemistry identifier.
        config: Source chemistry definition.
        project_root: Project root used to resolve reference paths.

    Returns:
        Compiled chemistry definition model.
    """
    # Compile libraries.
    libraries = {
        library_type: compile_library_info(library_type, library_config)
        for library_type, library_config in config[STRUCTURE].items()
    }

    available_feature_types = {
        features
        for library in libraries.values()
        for features in library.get_feature_types()
    }
    # Compile references.
    references = {
        feature_type: ReferenceInfo(
            feature_type=feature_type, path=Path(ref.get("path", None))
        )
        for feature_type, ref in (config.get(REFERENCES) or {}).items()
        if feature_type in available_feature_types
    }
    chemistry_def = ChemistryDef(
        id=chemistry_id,
        description=config.get(DESCRIPTION, "No description provided."),
        libraries=libraries,
        references=references,
    )
    validate_references(chemistry_def, project_root)
    demultiplexing_strategy = compile_demultiplexing_strategy(config, chemistry_def)
    chemistry_def.demultiplexing_strategy = demultiplexing_strategy
    return chemistry_def


def compile_library_info(
    library_type: str, library_config: dict[str, Any]
) -> LibraryInfo:
    """Compile one read structure into a runtime library model.

    Args:
        library_type: Source read structure key, for example `rna` or `dna`.
        library_config: Read structure configuration.

    Returns:
        Compiled library info.
    """
    pattern = library_config[PATTERN]
    sequence_parts = build_sequence_parts(library_config[MAPPINGS])
    validate_pattern(pattern, sequence_parts)
    update_sequence_parts_with_positions(
        sequence_parts=sequence_parts,
        pattern=pattern,
        offset=library_config.get(OFFSET, 0),
    )
    validated_feature_types = set()
    feature_types = {part.type for part in sequence_parts}
    for feature_type in feature_types:
        if feature_type not in LIST_OF_FEATURE_CLASSES:
            raise ValueError(
                f"Library type '{library_type}' has invalid feature type '{feature_type}'.\nDid you update the Feature class?"
            )
        else:
            validated_feature_types.add(feature_type)
    return LibraryInfo(
        library_type=library_type,
        sequence_parts=sorted(
            sequence_parts,
            key=lambda part: part.start if part.start is not None else -1,
        ),
    )


def build_sequence_parts(mappings: dict[str, dict[str, str]]) -> list[SequencePart]:
    """Build sequence parts from source feature mappings.

    Args:
        mappings: Mapping of feature type to symbol/name pairs.

    Returns:
        Sequence part models without positions.
    """
    sequence_parts: list[SequencePart] = []
    for part_type, part_info in mappings.items():
        for letter, name in part_info.items():
            sequence_parts.append(
                SequencePart(
                    letter=letter,
                    name=name,
                    type=part_type,
                )
            )
    return sequence_parts


def update_sequence_parts_with_positions(
    sequence_parts: list[SequencePart], pattern: str, offset: int
) -> None:
    """Add start and length coordinates to sequence parts.

    Args:
        sequence_parts: Sequence parts to update in place.
        pattern: Symbolic read pattern.
        offset: Number of bases to skip before parsing.
    """
    lengths = Counter(pattern)
    for sequence_part in sequence_parts:
        sequence_part.length = lengths[sequence_part.letter]
        sequence_part.start = pattern.index(sequence_part.letter) + offset


def validate_pattern(pattern: str, sequence_parts: list[SequencePart]) -> None:
    """Validate that a pattern and mappings describe the same symbols.

    Args:
        pattern: Symbolic read pattern.
        sequence_parts: Parts built from mappings.

    Raises:
        ValueError: If mappings and pattern disagree.
    """
    pattern_symbols = set(pattern)
    mapped_symbols = {sequence_part.letter for sequence_part in sequence_parts}
    missing_symbols = pattern_symbols - mapped_symbols
    if missing_symbols:
        raise ValueError(
            "Mapping keys do not cover all structure characters. "
            f"Missing: {','.join(sorted(missing_symbols))}"
        )
    extra_symbols = mapped_symbols - pattern_symbols
    if extra_symbols:
        raise ValueError(
            "Mapping keys contain extra characters not in structure. "
            f"Extra: {','.join(sorted(extra_symbols))}"
        )
    validate_unique_stretches(pattern)


def validate_unique_stretches(pattern: str) -> None:
    """Validate that each pattern symbol appears in exactly one contiguous stretch."""
    seen: set[str] = set()
    previous = ""
    for char in pattern:
        if char != previous:
            if char in seen:
                raise ValueError(
                    f"Character '{char}' appears in multiple stretches in pattern '{pattern}'."
                )
            seen.add(char)
            previous = char


def validate_references(chemistry_def: ChemistryDef, project_root: Path) -> None:
    """Validate declared reference files against mapped feature symbols. Check for file size and required columns.

    Args:
        chemistry_def: Compiled chemistry model.
        project_root: Project root used to resolve relative paths.
    """
    for feature_type, reference in chemistry_def.references.items():
        reference_path = resolve_reference_path(project_root, reference.path)
        if reference_path is None:
            continue
        ref_size = reference_path.stat().st_size
        if ref_size > 100 * 1024 and reference_path.suffix not in [".parquet"]:
            raise ValueError(
                f"Reference '{reference.path}' for chemistry '{chemistry_def.id}' is larger than 100kb ({ref_size} bytes). Convert file to .parquet"
            )
        reference_df = read_reference(reference_path)
        required_columns = REQUIRED_COLUMNS_BY_FEATURE_TYPE.get(feature_type, set())
        missing_columns = required_columns - set(reference_df.columns)
        if missing_columns:
            raise ValueError(
                f"Reference '{reference.path}' for chemistry '{chemistry_def.id}' is "
                f"missing required columns: {sorted(missing_columns)}"
            )

        expected_symbols = get_expected_symbols(chemistry_def, feature_type)
        if reference_df.height == 0:
            continue
        available_symbols = set(reference_df[SYMBOL].drop_nulls().to_list())
        missing_symbols = expected_symbols - available_symbols
        if missing_symbols:
            raise ValueError(
                f"Reference '{reference.path}' for chemistry '{chemistry_def.id}' is "
                f"missing symbols: {sorted(missing_symbols)}"
            )


def validate_demultiplexing_symbols(
    chemistry_def: ChemistryDef, demultiplexing_parts: list[str] | None
) -> None:
    """Validate demultiplexing symbols against barcode mappings."""
    if not demultiplexing_parts:
        return
    available_feature_types = chemistry_def.get_available_feature_types()
    for feature_type in available_feature_types:
        symbols = get_expected_symbols(chemistry_def, feature_type)
        if symbols.intersection(demultiplexing_parts):
            return

    raise ValueError(
        f"Chemistry '{chemistry_def.id}' is missing the demultiplexing symbols '{demultiplexing_parts}' in its mappings."
    )


def compile_demultiplexing_strategy(
    config, chemistry_def: ChemistryDef
) -> DemultiplexingStrategy | None:
    demultiplexing_parts = get_demultiplexing_parts_from_config(config)
    validate_demultiplexing_symbols(
        chemistry_def,
        demultiplexing_parts=demultiplexing_parts,
    )
    demultiplexing_in_n = symbol_in_libraries(chemistry_def, demultiplexing_parts)
    demultiplexing_strategy = deduct_demultiplexing_type(
        chemistry_def, demultiplexing_in_n, demultiplexing_parts
    )
    return demultiplexing_strategy


def get_demultiplexing_parts_from_config(config: dict[str, Any]) -> list[str]:
    """Return demultiplexing symbols from source config."""
    demultiplexing = config.get(DEMULTIPLEXING) or {}
    return demultiplexing.get(PARTS) or []


def get_expected_symbols(chemistry_def: ChemistryDef, feature_type: str) -> set[str]:
    """Return all symbols mapped to a feature type across libraries."""
    return {
        part.letter
        for library in chemistry_def.libraries.values()
        for part in library.sequence_parts
        if part.type.lower() == feature_type.lower()
    }


def symbol_in_libraries(chemistry_def: ChemistryDef, symbol: list[str]) -> int:
    """Check in how many libraries a symbol is mapped"""
    count = sum(
        1
        for library in chemistry_def.libraries.values()
        for part in library.sequence_parts
        if part.letter in symbol
    )
    return count


def resolve_reference_path(
    project_root: Path, reference_path: Path | None
) -> Path | None:
    """Resolve a reference path relative to the project root."""
    if reference_path is None:
        return None
    if reference_path.is_absolute():
        return reference_path
    return project_root / reference_path


def read_reference(reference_path: Path) -> pl.DataFrame:
    """Read a chemistry reference file.

    Args:
        reference_path: CSV, CSV.GZ, or Parquet reference path.

    Raises:
        ValueError: If the file is missing or has an unsupported extension.

    Returns:
        Reference data frame.
    """
    if not reference_path.exists():
        raise ValueError(f"Reference path {reference_path} does not exist.")

    suffixes = reference_path.suffixes
    extension = ".gz" if suffixes[-2:] == [".csv", ".gz"] else reference_path.suffix
    reader = REFERENCE_READERS.get(extension)
    if reader == "csv":
        return pl.read_csv(reference_path)
    if reader == "parquet":
        return pl.read_parquet(reference_path)
    raise ValueError(
        f"Unsupported reference format for {reference_path}. "
        "Expected .csv, .csv.gz, or .parquet."
    )
