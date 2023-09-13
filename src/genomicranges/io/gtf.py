import logging
from typing import Dict

from pandas import DataFrame, read_csv

# from ..GenomicRanges import GenomicRanges
from .pdf import from_pandas

# Variation of https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _parse_all_attribute(row: str) -> Dict:
    """Extract all keys from the gtf/gff attribute string.

    Args:
        row (str): A row from GTF.

    Returns:
        Dict: A dictionary containing extracted keys and their values.
    """
    attr = row["group"]
    infos = attr.split(";")
    vals = {}
    for info in infos:
        if info == "" or len(info) == 0:
            continue
        imap = info.strip().split(" ", 1)
        vals[imap[0]] = imap[1].strip().strip('"')
    return {**row, **vals}


def parse_gtf(path: str, compressed: bool) -> DataFrame:
    """Read a GTF file as :py:class:`~pandas.DataFrame`.

    Args:
        path (str): Path to the GTF file.
        compressed (bool): Whether the file is gzip compressed.

    Returns:
        DataFrame: Genome annotations from GTF as pandas dataframe.
    """

    from joblib import Parallel, delayed

    logging.info(f"Reading File - {path}")
    if compressed:
        df = read_csv(
            path,
            sep="\t",
            names=[
                "seqnames",
                "source",
                "feature",
                "starts",
                "ends",
                "score",
                "strand",
                "frame",
                "group",
            ],
            compression="gzip",
        )
    else:
        df = read_csv(
            path,
            sep="\t",
            names=[
                "seqnames",
                "source",
                "feature",
                "starts",
                "ends",
                "score",
                "strand",
                "frame",
                "group",
            ],
        )

    rows = Parallel(n_jobs=-2)(
        delayed(_parse_all_attribute)(row) for _, row in df.iterrows()
    )
    gtf = DataFrame.from_records(rows)
    gtf.drop(["group"], axis=1)

    return gtf


def read_gtf(file: str) -> "GenomicRanges":
    """Read  GTF file as :py:class:`~genomicranges.GenomicRanges.GenomicRanges`.

    Args:
        file (str): Path to GTF file.

    Returns:
        GenomicRanges:  Genome annotations from GTF.
    """
    compressed = True if file.endswith("gz") else False
    data = parse_gtf(file, compressed=compressed)

    return from_pandas(data)
