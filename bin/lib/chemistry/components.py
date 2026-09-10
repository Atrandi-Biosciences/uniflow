from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
from pathlib import Path
from lib.chemistry.common_const import LETTER, START, LENGTH, TYPE


@dataclass
class SequencePart:
    """A named segment extracted from a sequencing read.

    Attributes:
        letter: Symbol used in the chemistry pattern.
        name: Human-readable output column name.
        type: Feature group, for example barcode, handle, linker, or UMI.
        start: Zero-based start coordinate after applying read offset.
        length: Segment length in bases.
    """

    letter: str
    name: str
    type: str
    start: int | None = None
    length: int | None = None

    def __post_init__(self) -> None:
        """Infer correction behavior when it is not provided explicitly."""

    def to_dict(self) -> dict[str, dict[str, Any]]:
        """Serialize the sequence part to the runtime JSON shape.

        Returns:
            A single-item mapping keyed by the output column name.
        """
        return {
            self.name: {
                LETTER: self.letter,
                TYPE: self.type,
                START: self.start,
                LENGTH: self.length,
            }
        }


@dataclass
class LibraryInfo:
    """Compiled structure for one read library type."""

    library_type: str
    sequence_parts: list[SequencePart]

    def to_dict(self) -> list[dict[str, dict[str, Any]]]:
        """Serialize the library to the runtime JSON shape.

        Returns:
            Ordered sequence part dictionaries.
        """
        return [part.to_dict() for part in self.sequence_parts]

    def get_ordered_letters_by_type(self, feature_type: str) -> list[str]:
        """Return letters for parts matching a feature type in extraction order.

        Args:
            feature_type: Feature group to match case-insensitively.

        Returns:
            Pattern letters for matching parts.
        """
        return [
            part.letter for part in self.sequence_parts if part.type == feature_type
        ]

    def get_feature_types(self) -> set[str]:
        """Return the set of feature types used in this library."""
        return {part.type for part in self.sequence_parts}


@dataclass(frozen=True)
class ReferenceInfo:
    """External sequence reference declared by a chemistry."""

    feature_type: str
    path: Path | None

    def to_dict(self) -> dict[str, str]:
        """Serialize the reference path.

        Returns:
            A dictionary containing the source path as declared in YAML.
        """
        if self.path is None:
            return {}
        return {"path": self.path.as_posix()}


@dataclass
class DemultiplexingStrategy:
    """Inferred strategy for demultiplexing based on barcode configuration."""

    type: str | None = None
    symbols: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Serialize the demultiplexing strategy for runtime consumers."""
        return {TYPE: self.type, "symbols": self.symbols}
