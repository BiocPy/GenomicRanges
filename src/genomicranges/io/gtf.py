import logging
from typing import Dict, List, Union

# Variation of https://github.com/epiviz/epivizfileserver/src/epivizfileserver/cli.py

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _parse_all_attribute(row: str) -> Dict:
    """Extract all keys from the gtf/gff attribute string.

    Args:
        row:
            A row from GTF.

    Returns:
        A dictionary containing extracted keys and their values.
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


def parse_gtf(
    path: str,
    compressed: bool,
    skiprows: Union[int, List[int]] = None,
    comment: str = "#",
):
    """Read a GTF file as :py:class:`~pandas.DataFrame`.

    Args:
        path:
            Path to the GTF file.

        compressed:
            Whether the file is gzip compressed.

        skiprows:
            Rows to skip if the gtf file has header.

        comment:
            Character indicating that the line should not be
            parsed. Defaults to "#".

    Returns:
        Pandas DataFrame containing annotations from GTF.
    """

    from joblib import Parallel, delayed
    from pandas import DataFrame, read_csv

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
            skiprows=skiprows,
            comment=comment,
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
            skiprows=skiprows,
            comment=comment,
        )

    rows = Parallel(n_jobs=-2)(
        delayed(_parse_all_attribute)(row) for _, row in df.iterrows()
    )
    gtf = DataFrame.from_records(rows)
    gtf.drop(["group"], axis=1)

    return gtf


def read_gtf(
    file: str,
    skiprows: Union[int, List[int]] = None,
    comment: str = "#",
) -> "GenomicRanges":
    """Read a GTF file as :py:class:`~genomicranges.GenomicRanges.GenomicRanges`.

    Args:
        file:
            Path to GTF file.

        skiprows:
            Rows to skip if the gtf file has header.

        comment:
            Character indicating that the line should not be
            parsed. Defaults to "#".

    Returns:
        Genomic Ranges with annotations from the GTF file.
    """
    compressed = True if file.endswith("gz") else False
    data = parse_gtf(file, compressed=compressed, skiprows=skiprows, comment=comment)
    from ..GenomicRanges import GenomicRanges

    return GenomicRanges.from_pandas(data)
