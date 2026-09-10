# lib/pipeline/test_record_new.py

from __future__ import annotations

import __main__

import pytest
from pydantic import ValidationError

from lib.common_const import LibraryType, Modality
from lib.pipeline.record import PipelineRecord


@pytest.fixture
def library_type() -> LibraryType:
    """Return a valid LibraryType without depending on member names."""
    return next(iter(LibraryType))


@pytest.fixture
def modality() -> Modality:
    """Return a valid Modality without depending on member names."""
    return next(iter(Modality))


@pytest.fixture
def record(
    library_type: LibraryType,
    modality: Modality,
) -> PipelineRecord:
    return PipelineRecord(
        source_id="sample-001",
        library_type=library_type,
        modality=modality,
        process_name="unit-test",
    )


def test_initializes_with_expected_values(
    library_type: LibraryType,
    modality: Modality,
) -> None:
    record = PipelineRecord(
        source_id="sample-001",
        library_type=library_type,
        modality=modality,
        process_name="alignment",
    )

    assert record.source_id == "sample-001"
    assert record.library_type == library_type
    assert record.modality == modality
    assert record.process_name == "alignment"


def test_default_process_name_comes_from_main_file(
    monkeypatch: pytest.MonkeyPatch,
    library_type: LibraryType,
    modality: Modality,
) -> None:
    monkeypatch.setattr(
        __main__,
        "__file__",
        "/some/path/process_sample.py",
        raising=False,
    )

    record = PipelineRecord(
        source_id="sample-001",
        library_type=library_type,
        modality=modality,
    )

    assert record.process_name == "process_sample"


def test_default_process_name_is_interactive_when_main_file_is_missing(
    monkeypatch: pytest.MonkeyPatch,
    library_type: LibraryType,
    modality: Modality,
) -> None:
    monkeypatch.delattr(__main__, "__file__", raising=False)

    record = PipelineRecord(
        source_id="sample-001",
        library_type=library_type,
        modality=modality,
    )

    assert record.process_name == "interactive"


def test_rejects_extra_fields(
    library_type: LibraryType,
    modality: Modality,
) -> None:
    with pytest.raises(ValidationError, match="extra"):
        PipelineRecord(
            source_id="sample-001",
            library_type=library_type,
            modality=modality,
            unexpected_field="value",
        )


def test_validates_assignment(record: PipelineRecord) -> None:
    with pytest.raises(ValidationError):
        record.library_type = "not-a-library-type"
