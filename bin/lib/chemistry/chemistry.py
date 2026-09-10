from dataclasses import dataclass, field
from pathlib import Path
import polars as pl
from lib.chemistry.common_const import DESCRIPTION, START, LENGTH, LETTER, TYPE
from lib.common_const import SEQUENCE, READ_NAME
from typing import Any


from lib.chemistry.components import (
    DemultiplexingStrategy,
    LibraryInfo,
    ReferenceInfo,
    SequencePart,
)

from lib.chemistry.models import (
    FeatureBarcode,
)


@dataclass
class ChemistryDef:
    """Compiled chemistry definition used by the pipeline."""

    id: str
    description: str
    libraries: dict[str, LibraryInfo] = field(default_factory=dict)
    references: dict[str, ReferenceInfo] = field(default_factory=dict)
    demultiplexing_strategy: DemultiplexingStrategy | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the compiled chemistry to JSON-compatible data.

        Returns:
            Runtime-compatible chemistry configuration.
        """
        return {
            "id": self.id,
            "description": self.description,
            "libraries": {
                library_type: library_info.to_dict()
                for library_type, library_info in self.libraries.items()
            },
            "references": {
                feature_type: reference.to_dict()
                for feature_type, reference in self.references.items()
            },
            "demultiplexing_strategy": (
                self.demultiplexing_strategy.to_dict()
                if self.demultiplexing_strategy is not None
                else None
            ),
        }

    def load_from_config(self, config: dict[str, Any]) -> None:
        """Populate this model from a compiled runtime configuration.

        Args:
            config: One compiled chemistry configuration.
        """
        self.id = config.get("id", self.id)
        self.description = config.get(DESCRIPTION, self.description)
        self.libraries = {}
        for library_type, sequence_parts_config in (
            config.get("libraries", {}) or {}
        ).items():
            sequence_parts: list[SequencePart] = []
            for part_config in sequence_parts_config:
                for name, values in part_config.items():
                    sequence_parts.append(
                        SequencePart(
                            name=name,
                            letter=values.get(LETTER, ""),
                            type=values.get(TYPE, ""),
                            start=values.get(START, 0),
                            length=values.get(LENGTH, 0),
                        )
                    )
            self.libraries[library_type] = LibraryInfo(library_type, sequence_parts)

        self.references = {
            feature_type: ReferenceInfo(feature_type, Path(reference["path"]))
            for feature_type, reference in (config.get("references", {}) or {}).items()
        }

    @classmethod
    def from_general_config(
        cls, chemistry_id: str, general_config: dict[str, Any]
    ) -> "ChemistryDef":
        """Create a compiled chemistry model from a top-level config mapping.

        Args:
            chemistry_id: Chemistry identifier to load.
            general_config: Mapping of chemistry ids to compiled configs.

        Returns:
            Loaded chemistry definition.
        """
        chemistry_config = general_config.get(chemistry_id, {})
        instance = cls(
            id=chemistry_config.get("id", chemistry_id),
            description=chemistry_config.get("description", ""),
        )
        instance.load_from_config(chemistry_config)
        return instance

    def get_barcode_names(self, library_type: str) -> list[str]:
        """Return barcode output column names for a library type."""
        if library_type not in self.libraries:
            return []
        return [
            part.name
            for part in self.libraries[library_type].sequence_parts
            if part.type == FeatureBarcode.id
        ]

    def get_name_letter_mappings(self) -> dict[str, str]:
        """Return a mapping from output column names to pattern letters."""
        return {
            part.name: part.letter
            for library_info in self.libraries.values()
            for part in library_info.sequence_parts
        }

    def get_letter_name_mappings(self) -> dict[str, str]:
        """Return a mapping from pattern letters to output column names."""
        return {
            part.letter: part.name
            for library_info in self.libraries.values()
            for part in library_info.sequence_parts
        }

    def build_extract_exprs(self, library_type: str) -> list[tuple[str, pl.Expr]]:
        """Build Polars expressions to extract configured sequence parts.

        Args:
            library_type: Library type to extract.

        Raises:
            ValueError: If the library is unknown or a part lacks coordinates.

        Returns:
            Tuples of output column name and Polars expression.
        """
        if library_type not in self.libraries:
            raise ValueError(f"Unknown library type '{library_type}'")

        exprs: list[tuple[str, pl.Expr]] = []
        for part in self.libraries[library_type].sequence_parts:
            if part.start is None or part.length is None:
                raise ValueError(
                    f"Part '{part.name}' is missing 'start' or 'length' information"
                )
            exprs.append(
                (
                    part.name,
                    pl.col(SEQUENCE)
                    .str.slice(part.start, part.length)
                    .str.to_uppercase()
                    .alias(part.name),
                )
            )
        return exprs

    def get_available_feature_types(self) -> set[str]:
        """Return the set of all feature types used in this chemistry."""
        feature_types = set()
        for library_info in self.libraries.values():
            for part in library_info.sequence_parts:
                feature_types.add(part.type)
        return feature_types

    def build_extract_pipeline(
        self, lf: pl.LazyFrame, library_type: str
    ) -> pl.LazyFrame:
        """Build a LazyFrame pipeline that extracts sequence parts.

        Args:
            lf: LazyFrame containing read name and sequence columns.
            library_type: Library type to extract.

        Returns:
            LazyFrame with read name, sequence, and extracted columns.
        """
        select_exprs: list[pl.Expr] = [pl.col(READ_NAME), pl.col(SEQUENCE)]
        select_exprs.extend(expr for _, expr in self.build_extract_exprs(library_type))
        return lf.select(select_exprs)
