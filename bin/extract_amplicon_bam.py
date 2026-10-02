#!/usr/bin/env python3

import sys
import polars as pl
import polars_bio as plb


from lib.common_const import LibraryType, CHROM, READ_NAME

bam_path = sys.argv[1]
sample_name = sys.argv[2]
threads = int(sys.argv[3])


bam_plb_lf: pl.LazyFrame = (
    plb.scan_bam(
        str(bam_path), thread_num=threads, chunk_size=64, concurrent_fetches=threads
    )
    .with_columns(pl.col(CHROM).str.replace("^chr", "").str.to_uppercase())
    .rename({"name": READ_NAME})
    .drop(["mate_start", "mate_chrom"], strict=False)
)

bam_plb_lf.sink_parquet(f"{LibraryType.DNA.value}_bam.parquet")
