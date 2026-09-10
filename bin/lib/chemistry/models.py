"""Runtime models for compiled chemistry definitions."""

from lib.common_const import SEQUENCE
from lib.chemistry.common_const import WELL_COORDINATE, SYMBOL

from typing import ClassVar


class Feature:
    """Base class for features in a library. Each feature type should inherit from this class and have a unique 'id' attribute."""

    id: ClassVar[str] = ""
    required_columns: ClassVar[set[str]] = set()


class FeatureBarcode(Feature):
    """Feature representing a barcode in a library."""

    id: ClassVar[str] = "barcode"
    required_columns: ClassVar[set[str]] = {SYMBOL, SEQUENCE, WELL_COORDINATE}


class FeatureHandle(Feature):
    """Feature representing a handle in a library."""

    id: ClassVar[str] = "handle"
    required_columns: ClassVar[set[str]] = {SYMBOL, SEQUENCE}


class FeatureLinker(Feature):
    """Feature representing a linker in a library."""

    id: ClassVar[str] = "linker"
    required_columns: ClassVar[set[str]] = {SYMBOL, SEQUENCE}


class FeatureUMI(Feature):
    """Feature representing a UMI in a library."""

    id: ClassVar[str] = "UMI"
    required_columns: ClassVar[set[str]] = set()


class FeatureSampleIndex(Feature):
    """Feature representing a sample index in a library."""

    id: ClassVar[str] = "sample_index"
    required_columns: ClassVar[set[str]] = {SYMBOL, SEQUENCE}


LIST_OF_FEATURE_CLASSES: list[str] = [
    feature_type.id for feature_type in Feature.__subclasses__()
]

REQUIRED_COLUMNS_BY_FEATURE_TYPE: dict[str, set[str]] = {
    feature_type.id: feature_type.required_columns
    for feature_type in Feature.__subclasses__()
}
