from typing import Final

START: Final[str] = "start"
LENGTH: Final[str] = "length"
LETTER: Final[str] = "letter"
SYMBOL: Final[str] = "symbol"

STRUCTURE: Final[str] = "structure"
REFERENCES: Final[str] = "references"
DEMULTIPLEXING: Final[str] = "demultiplexing"
MAPPINGS: Final[str] = "mappings"
OFFSET: Final[str] = "offset"
PATTERN: Final[str] = "pattern"
PARTS: Final[str] = "parts"
DESCRIPTION: Final[str] = "description"
WELL_COORDINATE: Final[str] = "well_coordinate"
TYPE: Final[str] = "type"

REFERENCE_READERS: Final[dict[str, str]] = {
    ".csv": "csv",
    ".gz": "csv",
    ".parquet": "parquet",
}
