#!/usr/bin/env python3
import sys
from pathlib import Path

import yaml
from lib.metrics.report import build_config

metrics_reference_csv: Path = Path(sys.argv[1])
id: str = sys.argv[2]
aggregated_metrics_csv: Path = Path(sys.argv[3])
scope: str = sys.argv[4]
out_config: Path = Path(sys.argv[5])

config = build_config(metrics_reference_csv, aggregated_metrics_csv, scope)

out_config.write_text(yaml.safe_dump(config, sort_keys=False))
print(f"Wrote {out_config}")
