#!/usr/bin/env python3
from lib.bam_tools.utils import extract_bam_parallel
import sys
import pysam
from pathlib import Path
import pyarrow as pa

# constants
from lib.bam_tools.utils import (
    GENE_SYMBOL_TAG,
    CORRECTED_UMI_TAG,
    CORRECTED_BARCODE_TAG,
)

bam_path: Path = Path(sys.argv[1])
sample_name: str = sys.argv[2]
library_type: str = sys.argv[3]
n_threads: int = int(sys.argv[4])
bam_index_path = bam_path.with_suffix(".bam.bai")
if not bam_index_path.exists():
    print(f"Indexing {bam_path} file")
    pysam.index(str(bam_path))

tags = [
    pa.field(GENE_SYMBOL_TAG, pa.string()),
    pa.field(CORRECTED_UMI_TAG, pa.string()),
    pa.field(CORRECTED_BARCODE_TAG, pa.string()),
]

extract_bam_parallel(
    bam_path=str(bam_path),
    out_parquet=f"{library_type}_bam.parquet",
    tags=tags,
    workers=n_threads,
)
