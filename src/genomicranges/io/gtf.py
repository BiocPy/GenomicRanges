from typing import MutableMapping
import pandas as pd
import logging
from joblib import Parallel, delayed

from ..GenomicRanges import GenomicRanges
from .pdf import fromPandas

## Variation of https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def parse_all_attribute(row: str) -> MutableMapping:
    """Extract all keys from the gtf/gff attribute string

    Args:
        row (str): a row from GTF

    Returns:
        MutableMapping: dict containing extracted keys and their values
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


def parse_gtf(path: str, compressed: bool) -> pd.DataFrame:
    """Read a GTF file as pandas DataFrame

    Args:
        path (str): path to the file
        compressed (bool): is compression gzip ?

    Returns:
        pd.DataFrame: DataFrame representation of the gtf file
    """
    logging.info(f"Reading File - {path}")
    if compressed:
        df = pd.read_csv(
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
        df = pd.read_csv(
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
        delayed(parse_all_attribute)(row) for _, row in df.iterrows()
    )
    gtf = pd.DataFrame.from_records(rows)
    gtf.drop(["group"], axis=1)

    return gtf


def readGTF(file: str) -> GenomicRanges:
    """Read  GTF file as `GenomicRanges`.

    Args:
        file (str): path to gtf file

    Returns:
        GenomicRanges:  a new `GenomicRanges` with the genomic regions.
    """
    compressed = True if file.endswith("gz") else False
    data = parse_gtf(file, compressed=compressed)

    return fromPandas(data)
