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


def test_reduce():
    assert gr is not None

    reduced_gr = gr.reduce()

    assert reduced_gr is not None
    assert reduced_gr.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (reduced_gr.start == np.array([101, 104, 102, 105, 103])).all()
    assert (reduced_gr.width == np.array([11, 30, 21, 5, 25])).all()
    assert (reduced_gr.strand == np.array([0, 1, -1, -1, 0])).all()

    reduced_gr = gr.reduce(ignore_strand=True)

    assert reduced_gr is not None
    assert reduced_gr.seqnames == ["chr1", "chr2", "chr3"]
    assert (reduced_gr.start == np.array([101, 102, 103])).all()
    assert (reduced_gr.width == np.array([11, 32, 25])).all()
    assert (reduced_gr.strand == np.array([0, 0, 0])).all()


def test_reduce_with_gapwidth():
    assert gr is not None

    reduced_gr = gr.reduce(min_gap_width=10)

    assert reduced_gr is not None
    assert reduced_gr.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (reduced_gr.start == np.array([101, 104, 102, 105, 103])).all()
    assert (reduced_gr.width == np.array([11, 30, 21, 5, 25])).all()
    assert (reduced_gr.strand == np.array([0, 1, -1, -1, 0])).all()


def test_reduce_with_gapwidth_with_reverse_map():
    assert gr is not None

    reduced_gr = gr.reduce(min_gap_width=10, with_reverse_map=True)

    assert reduced_gr is not None
    assert reduced_gr.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (reduced_gr.start == np.array([101, 104, 102, 105, 103])).all()
    assert (reduced_gr.width == np.array([11, 30, 21, 5, 25])).all()
    assert (reduced_gr.strand == np.array([0, 1, -1, -1, 0])).all()
    assert reduced_gr.mcols.column("revmap") == [[0], [3], [1], [4], [2]]


def test_range():
    assert gr is not None

    range_gr = gr.range()

    assert range_gr is not None
    assert range_gr.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (range_gr.start == np.array([101, 104, 102, 105, 103])).all()
    assert (range_gr.width == np.array([11, 30, 21, 5, 25])).all()
    assert (range_gr.strand == np.array([0, 1, -1, -1, 0])).all()

    range_gr = gr.range(ignore_strand=True)

    assert range_gr is not None
    assert range_gr.seqnames == ["chr1", "chr2", "chr3"]
    assert (range_gr.start == np.array([101, 102, 103])).all()
    assert (range_gr.width == np.array([11, 32, 25])).all()
    assert (range_gr.strand == np.array([0, 0, 0])).all()


def test_gaps():
    assert gr is not None

    print("gr seqnames", gr.seqnames, gr._seqnames)
    print("gr., seqinfo", gr.seqinfo.seqnames)

    out = gr.gaps()

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([1] * 5)).all()
    assert (out.width == np.array([100, 103, 101, 104, 102])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()


# def test_gaps_with_start():
#     assert gr is not None

#     out = gr.gaps(start=5)

#     assert out is not None
#     assert out.seqnames == ["chr1", "chr2", "chr2", "chr3"]
#     assert out.start == [5, 5, 104, 5]
#     assert out.width == [100, 101, 108, 103]
#     assert out.strand == ["*", "-", "-", "+"]


# def test_gaps_with_start_filter():
#     assert gr is not None

#     out = gr.gaps(start=103)

#     assert out is not None
#     assert out.seqnames == ["chr2", "chr3"]
#     assert out.start == [104, 103]
#     assert out.width == [108, 103]
#     assert out.strand == ["-", "+"]


# def test_gaps_with_end():
#     assert gr is not None

#     out = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

#     assert out is not None
#     assert out.seqnames == [
#         "chr1",
#         "chr1",
#         "chr1",
#         "chr2",
#         "chr2",
#         "chr2",
#         "chr2",
#         "chr2",
#         "chr3",
#         "chr3",
#         "chr3",
#     ]
#     assert out.start == [1, 1, 1, 1, 1, 1, 104, 112, 1, 1, 1]
#     assert out.width == [
#         100,
#         120,
#         120,
#         120,
#         120,
#         101,
#         108,
#         120,
#         120,
#         103,
#         120,
#     ]

#     assert out.strand == [
#         "*",
#         "+",
#         "-",
#         "*",
#         "+",
#         "-",
#         "-",
#         "-",
#         "*",
#         "+",
#         "-",
#     ]


# def test_disjoin():
#     assert gr is not None

#     disjoin_gr = gr.disjoin()

#     assert disjoin_gr is not None
#     assert disjoin_gr.column("seqnames") == [
#         "chr1",
#         "chr1",
#         "chr1",
#         "chr2",
#         "chr2",
#         "chr3",
#     ]
#     assert disjoin_gr.column("starts") == [101, 103, 113, 102, 109, 104]
#     assert disjoin_gr.column("ends") == [102, 112, 128, 103, 111, 134]
#     assert disjoin_gr.column("strand") == ["*", "*", "*", "-", "-", "+"]
