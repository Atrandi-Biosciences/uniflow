from __future__ import annotations
import __main__
from typing import Any
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict, Field

from lib.common_const import (
    LibraryType,
    Modality,
    SOURCE_ID,
    LIBRARY_TYPE,
    MODALITY,
    PROCESS_NAME,
)


def _default_process_name() -> str:
    """Build a default process name from the current Python entrypoint.

    Returns:
        str: The stem of ``__main__.__file__`` when running from a script,
        or ``"interactive"`` when no entrypoint file is available.
    """
    main_file = getattr(__main__, "__file__", None)

    if main_file is None:
        return "interactive"

    return Path(main_file).stem


class PipelineRecord(BaseModel):
    """Container for pipeline key/value records scoped to one sample context.

    The model stores record rows in an internal Polars DataFrame with a fixed
    schema (key and value columns). Metadata fields identifying the sample and
    processing context are added when exporting via :meth:`finalize_df`.
    """

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        arbitrary_types_allowed=True,
    )

    source_id: str
    library_type: LibraryType
    modality: Modality
    process_name: str = Field(default_factory=_default_process_name)
    _key_name: str = "key"
    _value_name: str = "value"
    _value_type = pl.String
    records_df: pl.DataFrame = Field(
        default_factory=pl.DataFrame,
        exclude=True,
    )

    def model_post_init(self, __context: Any) -> None:
        """Initialize the internal records DataFrame after model creation.

        Args:
            __context: Pydantic post-init context. Unused here.
        """
        self.records_df = pl.DataFrame(
            schema={
                self._key_name: pl.String,
                self._value_name: self._value_type,
            }
        )

    def add(self, key: str, value: Any) -> None:
        """Append a single record row to the in-memory table.

        Args:
            key: Record key name.
            value: Record value associated with ``key``.
        """
        key, value = self._validate_record(key, value)
        self.records_df = self.records_df.vstack(
            pl.DataFrame(
                {self._key_name: [key], self._value_name: [value]},
                schema={self._key_name: pl.String, self._value_name: self._value_type},
            )
        )

    def extend(self, records: pl.DataFrame) -> None:
        """Append multiple record rows from another DataFrame.

        Args:
            records: Polars DataFrame containing rows matching the internal
                key/value schema.

        Raises:
            TypeError: If ``records`` is not a Polars DataFrame.
            ValueError: If ``records`` columns do not match the internal
                records schema.
        """
        if not isinstance(records, pl.DataFrame):
            raise TypeError(
                f"Expected a Polars DataFrame, got {type(records).__name__}"
            )
        if records.columns != self.records_df.columns:
            raise ValueError(
                f"DataFrame columns {records.columns} do not match expected columns {self.records_df.columns}"
            )
        self.records_df = self.records_df.vstack(records)

    def _validate_record(self, key: str, value: Any) -> tuple[str, Any]:
        """Validate or normalize a record before insertion.

        Subclasses can override this hook to enforce allowed keys, cast values,
        or perform other record-level validation.

        Args:
            key: Proposed record key.
            value: Proposed record value.

        Returns:
            tuple[str, Any]: The validated ``(key, value)`` pair.
        """
        return key, value

    def write_records(self, path: Path | str, extension="parquet") -> None:
        """Write finalized records to disk as parquet or CSV.

        Args:
            path: Output path for the file.
            extension: Output format. Supported values are ``"parquet"`` and
                ``"csv"``.

        Raises:
            ValueError: If there are no records to write.
            ValueError: If ``extension`` is not supported.
        """
        if self.records_df.shape[0] == 0:
            return
        self.finalize_df()
        if extension == "parquet":
            self.records_df.write_parquet(path)
        elif extension == "csv":
            self.records_df.write_csv(path)
        else:
            raise ValueError(f"Unsupported extension: {extension}")

    def finalize_df(self) -> None:
        """Attach metadata columns and reorder fields for downstream output.

        This method mutates ``records_df`` in place, adding source and process
        metadata columns and selecting the canonical output column order.

        Raises:
            ValueError: If no records are present.
        """

        if self.records_df.shape[0] == 0:
            raise ValueError("No records to extract")
        self.records_df = (
            self.records_df.with_columns(
                pl.lit(self.source_id).alias(SOURCE_ID),
                pl.lit(self.library_type.value).alias(LIBRARY_TYPE),
                pl.lit(self.modality.value).alias(MODALITY),
                pl.lit(self.process_name).alias(PROCESS_NAME),
            )
            .select(
                SOURCE_ID,
                LIBRARY_TYPE,
                MODALITY,
                PROCESS_NAME,
                self._key_name,
                self._value_name,
            )
            .sort([SOURCE_ID, LIBRARY_TYPE, MODALITY, PROCESS_NAME, self._key_name])
        )
