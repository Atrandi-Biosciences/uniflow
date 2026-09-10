import pysam
from typing import List, Tuple, Optional, Final
from pathlib import Path
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import multiprocessing as mp
import os

from lib.common_const import MAPPING_QUALITY, READ_NAME

GENE_SYMBOL_TAG: Final[str] = "GN"
RAW_UMI_TAG: Final[str] = "UR"
CORRECTED_UMI_TAG: Final[str] = "UB"
CORRECTED_BARCODE_TAG: Final[str] = "CB"
STAR_NO_GENE: Final[str] = "-"
STAR_BAD_UMI: Final[str] = "-"
STAR_MIN_MAPPING_QUALITY: Final[int] = 255


def make_shards(bam_path: str, window_mb: int = 10) -> List[Tuple[str, int, int]]:
    """Build genomic (contig, start, end) shards from BAM header."""
    shards = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for contig, length in zip(bam.header.references, bam.header.lengths):
            win = window_mb * 1_000_000
            for s in range(0, length, win):
                e = min(s + win, length)
                shards.append((contig, s, e))
    shards.append(("*", None, None))  # Add unmapped reads shard
    return shards


# TODO: Write this function as a general parser with a provided dict for tags and flags.
def process_shard(args):
    """Process bam shards for data extraction and stream it to parquet files

    Args:
        args (_type_): _description_

    Returns:
        str: path to sharded parquet file
    """
    (
        bam_path,
        out_dir,
        shard,
        tags,
        batch_size,
        compression,
        clevel,
    ) = args
    contig, start, end = shard
    if contig == "*":
        out_path = Path(out_dir) / "part_unmapped.parquet"
    else:
        out_path = Path(out_dir) / f"part_{contig}_{start}_{end}.parquet"
    # Extract the additional tags
    tag_keys = [i.name for i in tags]
    default_fields = [
        pa.field(MAPPING_QUALITY, pa.int32()),
        pa.field(READ_NAME, pa.string()),
    ]

    tags.extend(default_fields)

    schema = pa.schema(i for i in tags)

    writer = pq.ParquetWriter(
        out_path,
        schema,
        compression=compression,
        use_dictionary=True,
        compression_level=clevel,
    )
    batch = {i.name: [] for i in tags}

    count = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Note: fetch is inclusive-exclusive on coordinates

        for aln in bam.fetch(contig, start, end):
            mq = int(aln.mapping_quality)
            qname = aln.query_name
            for key in tag_keys:
                try:
                    value = aln.get_tag(key)
                except KeyError:
                    value = None
                batch[key].append(value)

            batch[MAPPING_QUALITY].append(mq)
            batch[READ_NAME].append(qname)
            count += 1

            if count % batch_size == 0:
                writer.write_table(pa.table(batch, schema=schema))
                for k in batch:
                    batch[k].clear()

        if batch[READ_NAME]:
            writer.write_table(pa.table(batch, schema=schema))

    writer.close()
    return str(out_path)


def merge_parquets(parts: List[str], out_path: Path) -> None:
    """Merge parquet files into one

    Args:
        parts (List[str]): Individual parquet files
        out_path (str): output path of the concatenated parquet file
    """
    # Fast concat of row groups
    datasets = [pq.ParquetFile(p) for p in parts]
    writer = None
    for pf in datasets:
        for i in range(pf.num_row_groups):
            rg = pf.read_row_group(i)
            if writer is None:
                writer = pq.ParquetWriter(
                    out_path,
                    rg.schema,
                    compression=pf.metadata.row_group(0).column(0).compression,
                )
            writer.write_table(rg)
    if writer is not None:
        writer.close()


def extract_bam_parallel(
    bam_path: str,
    out_parquet: str,
    tags: list,
    window_mb: int = 10,
    workers: Optional[int] = None,
    batch_size: int = 200_000,
    compression: str = "zstd",
    compression_level: int = 3,
    tmp_dir: Optional[str] = None,
):
    """Extract read information from a BAM file coming from STAR in parallel by sharding chromosome chunks

    Args:
        bam_path (str): Path to the bam file
        out_parquet (str): output of the final parquet df
        tags (dict): list of tags to extract
        window_mb (int, optional): Windows size for chunking. Defaults to 10.
        workers (Optional[int], optional): Number of max workers. Defaults to None.
        batch_size (int, optional): Batch size of writing out. Defaults to 200_000.
        compression (str, optional): compression type. Defaults to "zstd".
        compression_level (int, optional): Compressino level. Defaults to 3.
        tmp_dir (Optional[str], optional): Temp dir, will write to dir of out_parquet if None. Defaults to None.
    """
    workers = workers or max(1, mp.cpu_count() // 2)
    tmp_dir = tmp_dir or (Path(out_parquet).with_suffix("").as_posix() + "_parts")
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    shards = make_shards(bam_path, window_mb)
    tasks = [
        (
            bam_path,
            tmp_dir,
            shard,
            tags,
            batch_size,
            compression,
            compression_level,
        )
        for shard in shards
    ]

    with mp.Pool(processes=workers) as pool:
        parts = list(pool.imap_unordered(process_shard, tasks, chunksize=1))

    # Merge
    tmp_concatenated = Path(tmp_dir) / "concatenated.parquet"
    merge_parquets(parts, tmp_concatenated)
    pl.scan_parquet(tmp_concatenated).unique().sink_parquet(out_parquet)
    parts.append(tmp_concatenated.as_posix())
    # Cleanup
    for p in parts:
        try:
            os.remove(p)
        except OSError:
            pass
    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass
