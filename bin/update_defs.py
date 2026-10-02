#!/usr/bin/env python3
"""Compile source chemistry definitions into runtime JSON files."""

from __future__ import annotations


import logging
import sys
from pathlib import Path

from lib.chemistry.compiler import compile_chemistry_definitions
from lib.chemistry.loader import (
    load_chemistry_yaml,
    validate_chemistry_yaml,
    write_chemistry_json,
)

chemistry_def_path = Path(sys.argv[1])
ASSETS_ROOT_PATH = Path(sys.argv[2])
schema_path = Path(sys.argv[3])

chemistry_defs_output = Path("chemistry_defs.json")

logging.basicConfig(
    level=logging.INFO,
    format="{levelname}:{name}:{message}",
    style="{",
)


logging.info("Loading chemistry definitions from %s", chemistry_def_path)
source_config = load_chemistry_yaml(chemistry_def_path)
logging.info("Validating chemistry definitions with %s", schema_path)
validate_chemistry_yaml(source_config, schema_path)
compiled_config = compile_chemistry_definitions(
    chemistry_definitions=source_config,
    project_root=ASSETS_ROOT_PATH,
)

logging.info("Writing compiled chemistry definitions to %s", chemistry_defs_output)
write_chemistry_json(compiled_config, chemistry_defs_output)
