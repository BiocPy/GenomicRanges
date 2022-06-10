from typing import Union
import pandas as pd
import logging
from joblib import Parallel, delayed

## Variation of https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def parse_attribute(attr: str, key: str) -> Union[None, str]:
    """Extract a key from the gtf/gff attribute string

    Args:
        attr (str): attribute string
        key (str): key to extract

    Returns:
        Union[None, str]: value if key exists else `None`
    """
    if key in attr:
        tstr = attr.split(key, 1)
        tstrval = tstr[1].split(";", 1)
        return tstrval[0][1:]
    else:
        return None


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
                "seqname",
                "source",
                "feature",
                "start",
                "end",
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
                "seqname",
                "source",
                "feature",
                "start",
                "end",
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
