import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random
from iranges import IRanges
from biocframe import BiocFrame
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

gr = GenomicRanges(
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
    ranges=IRanges(range(101, 111), [10] * 10),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame({"score": range(0, 10), "GC": [random() for _ in range(10)]}),
)


def test_nearest():
    assert gr is not None

    q_gr = GenomicRanges(
        seqnames=["chr1", "chr2", "chr3"],
        ranges=IRanges([200, 105, 1190], [3, 1, 10]),
    )

    query_hits = gr.nearest(q_gr)

    assert query_hits is not None
    assert query_hits == [[4], [3], []]


def test_precede():
    assert gr is not None

    q_gr = GenomicRanges(
        seqnames=["chr1", "chr2", "chr3"],
        ranges=IRanges([200, 105, 1190], [3, 1, 10]),
    )

    query_hits = gr.precede(q_gr)

    assert query_hits is not None
    assert query_hits == [[], [], []]


def test_follow():
    assert gr is not None

    q_gr = GenomicRanges(
        seqnames=["chr1", "chr2", "chr3"],
        ranges=IRanges([200, 105, 1190], [3, 1, 10]),
    )

    query_hits = gr.follow(q_gr)

    assert query_hits is not None
    assert query_hits == [[4], [], []]


def test_distance():
    query = GenomicRanges(["A", "A", "A"], ranges=IRanges([1, 2, 10], [5, 7, 2]))
    subject = GenomicRanges(["A", "A", "A"], ranges=IRanges([6, 5, 13], [5, 6, 3]))

    distance = subject.distance(query)

    assert distance is not None
    assert isinstance(distance, np.ndarray)
    assert distance.tolist() == [0, 0, 1]
