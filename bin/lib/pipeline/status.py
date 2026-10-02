from __future__ import annotations

import sys
from enum import StrEnum
from typing import Any, NoReturn

import polars as pl

from lib.pipeline.record import PipelineRecord


class ReasonCode(StrEnum):
    INSUFFICIENT_READS = "INSUFFICIENT_READS"
    LOW_COVERAGE = "LOW_COVERAGE"
    LOW_BARCODES = "LOW_BARCODES"
    LOW_FEATURES = "LOW_FEATURES"
    NO_CELLS = "NO_CELLS"
    INCONSISTENT_MERGE = "INCONSISTENT_MERGE"


class StatusRecord(PipelineRecord):
    """A mapping of reason codes to human-readable messages.
    args:
        source_id: The source ID of the sample.
        library_type: The library type of the sample.
        modality: The modality of the sample.
    """

    _key_name: str = "reason_code"
    _value_name: str = "message"
    _value_type: type = pl.String

    def model_post_init(self, __context: Any) -> None:
        self.records_df = pl.DataFrame(
            schema={
                self._key_name: pl.String,
                self._value_name: self._value_type,
            }
        )

    def _validate_record(
        self,
        key: str,
        value: Any,
    ) -> tuple[str, str]:
        try:
            reason_code = ReasonCode(key)
        except ValueError as exc:
            allowed = ", ".join(code.value for code in ReasonCode)
            raise ValueError(
                f"Invalid reason code {key!r}. Allowed values: {allowed}"
            ) from exc

        if not isinstance(value, str):
            raise TypeError(
                f"Status record values must be strings, got {type(value).__name__}"
            )

        return reason_code.value, value

    def record_and_exit(self, reason: ReasonCode, message: str) -> NoReturn:
        """Record why this step produced no output, then exit 0.

        Not a failure: a sample with nothing to attribute is a valid outcome of
        count_variant*.

        Pairs the write with the exit because forgetting write_records leaves the
        reason in the task log only, where nothing collects it.
        """
        self.add(reason.value, message)
        self.write_records("status.parquet")
        print(message)
        sys.exit(0)
