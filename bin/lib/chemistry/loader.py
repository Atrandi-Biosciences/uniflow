"""File IO helpers for chemistry definitions."""

import json
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError


def load_chemistry_yaml(path: Path) -> dict[str, Any]:
    """Load source chemistry definitions from YAML.

    Args:
        path: Path to a chemistry definition YAML file.

    Raises:
        ValueError: If the YAML file is empty.

    Returns:
        Parsed chemistry definitions.
    """
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        raise ValueError(f"Chemistry definition file {path} is empty.")
    return data


def write_chemistry_json(config: dict[str, Any], path: Path) -> None:
    """Write compiled chemistry definitions as formatted JSON.

    Args:
        config: JSON-serializable chemistry configuration.
        path: Output path.
    """
    with path.open("w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=4)
        handle.write("\n")


def load_chemistry_json(path: Path) -> dict[str, Any]:
    """Load compiled chemistry definitions from JSON.

    Args:
        path: JSON path.

    Returns:
        Parsed compiled chemistry definitions.
    """
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_chemistry_yaml(config: dict[str, Any], schema_path: Path) -> None:
    """Validate source chemistry definitions against a JSON schema.

    Args:
        config: Parsed source chemistry definitions.
        schema_path: Path to the chemistry JSON schema.

    Raises:
        ValueError: If the source config does not match the schema.
    """
    schema = load_chemistry_json(schema_path)
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(config), key=lambda error: error.path)
    if not errors:
        return

    formatted_errors = "\n".join(format_schema_error(error) for error in errors)
    raise ValueError(
        f"Chemistry definition YAML failed schema validation:\n{formatted_errors}"
    )


def format_schema_error(error: ValidationError) -> str:
    """Format a JSON schema validation error for CLI output.

    Args:
        error: Validation error returned by jsonschema.

    Returns:
        Human-readable validation error.
    """
    path = ".".join(str(part) for part in error.absolute_path)
    location = path if path else "<root>"
    return f"- {location}: {error.message}"
