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
        "chr3",
        "chr2",
        "chr3",
    ],
    ranges=IRanges([x for x in range(101, 106)], [11, 21, 25, 30, 5]),
    strand=["*", "-", "*", "+", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 5),
            "GC": [random() for _ in range(5)],
        }
    ),
)


def test_matches():
    assert gr is not None

    query_hits = gr.match(gr[2:5])

    assert query_hits is not None
    assert query_hits.column("subject_hits") == [[2], [3], [4]]


# def test_duplicated():
#     assert gr is not None

#     query_hits = gr.duplicated()

#     assert query_hits is not None
#     assert query_hits == [
#         False,
#         False,
#         False,
#         False,
#         False,
#         False,
#         False,
#         False,
#         False,
#         False,
#     ]


# def test_is_unsorted():
#     assert gr is not None

#     result = gr.is_unsorted()

#     assert result is True


# def test_order():
#     assert gr is not None

#     order = gr.order()

#     assert order == [4, 5, 0, 3, 1, 2, 6, 7, 8, 9]


# def test_sort():
#     assert gr is not None
#     result = gr.sort()

#     assert result is not None
#     assert result.shape == gr.shape

#     result = gr.sort(decreasing=True)

#     assert result is not None
#     assert result.shape == gr.shape


# def test_rank():
#     assert gr is not None
#     result = gr.rank()

#     assert result is not None
#     assert result == [2, 4, 5, 3, 0, 1, 6, 7, 8, 9]
