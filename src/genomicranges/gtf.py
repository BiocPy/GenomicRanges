from typing import Union
import pandas as pd
import logging
from joblib import Parallel, delayed

## Variation of https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def parse_all_attribute(row: str) -> dict:
    """Extract all keys from the gtf/gff attribute string

    Args:
        row (str): a row from GTF

    Returns:
        dict: map of keys and its values
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
