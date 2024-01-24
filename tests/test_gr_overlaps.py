import pytest
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges
from biocframe import BiocFrame
from iranges import IRanges
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

subject = GenomicRanges(
    seqnames=[
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr1",
        "chr1",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
    ],
    ranges=IRanges(range(1, 11), [10] * 10),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    ),
)

query = GenomicRanges(
    seqnames=["chr2"],
    ranges=IRanges([4], [4]),
    strand=["+"],
)


def test_find_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query)

    assert res is not None
    assert isinstance(res, list)
    assert res == [[1, 2]]


def test_find_overlaps_query_type():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query, query_type="within")

    assert res is not None
    assert res == [[1, 2]]


def test_count_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.count_overlaps(query)

    assert res is not None
    assert isinstance(res, list)
    assert res == [2]


def test_subset_by_overlaps():
    query2 = GenomicRanges(
        seqnames=["chr2", "chr4"],
        ranges=IRanges(start=[4, 4], width=[2, 2]),
        strand=["+", "-"],
    )

    assert subject is not None
    assert query2 is not None

    res = subject.subset_by_overlaps(query2)

    assert res is not None
    assert isinstance(res, GenomicRanges)
    assert len(res) == 2
