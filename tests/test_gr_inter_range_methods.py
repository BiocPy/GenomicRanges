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

    out = gr.reduce()

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()

    out = gr.reduce(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr3"]
    assert (out.start == np.array([101, 102, 103])).all()
    assert (out.width == np.array([11, 32, 25])).all()
    assert (out.strand == np.array([0, 0, 0])).all()


def test_reduce_with_contigs():
    gr2 = GenomicRanges(
        seqnames=[
            "chr1_gl123",
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

    assert gr2 is not None

    out = gr2.reduce()

    assert out is not None
    assert out.seqnames == ["chr1_gl123", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()

    out = gr2.reduce(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1_gl123", "chr2", "chr3"]
    assert (out.start == np.array([101, 102, 103])).all()
    assert (out.width == np.array([11, 32, 25])).all()
    assert (out.strand == np.array([0, 0, 0])).all()


def test_reduce_with_gapwidth():
    assert gr is not None

    out = gr.reduce(min_gap_width=10)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()


def test_reduce_with_gapwidth_with_reverse_map():
    assert gr is not None

    out = gr.reduce(min_gap_width=10, with_reverse_map=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()
    assert out.mcols.column("revmap") == [[0], [3], [1], [4], [2]]


def test_range():
    assert gr is not None

    out = gr.range()

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()

    out = gr.range(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr3"]
    assert (out.start == np.array([101, 102, 103])).all()
    assert (out.width == np.array([11, 32, 25])).all()
    assert (out.strand == np.array([0, 0, 0])).all()


def test_gaps():
    assert gr is not None

    out = gr.gaps()
    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([1] * 5)).all()
    assert (out.width == np.array([100, 103, 101, 104, 102])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()

    out = gr.gaps(ignore_strand=True)
    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr3"]
    assert (out.start == np.array([1] * 3)).all()
    assert (out.width == np.array([100, 101, 102])).all()
    assert (out.strand == np.array([0, 0, 0])).all()


def test_gaps_with_start():
    assert gr is not None

    out = gr.gaps(start=5)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([5] * 5)).all()
    assert (out.width == np.array([96, 99, 97, 100, 98])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()


def test_gaps_with_start_filter():
    assert gr is not None

    out = gr.gaps(start=103)

    assert out is not None
    assert out.seqnames == ["chr2", "chr3"]
    assert (out.start == np.array([103] * 2)).all()
    assert (out.width == np.array([1, 2])).all()
    assert (out.strand == np.array([1, -1])).all()


def test_gaps_with_end():
    assert gr is not None

    out = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

    assert out is not None
    assert out.seqnames == [
        "chr1",
        "chr1",
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
    ]
    assert (out.start == np.array([1, 1, 1, 112, 1, 1, 1, 1, 1, 110, 1])).all()
    assert (
        out.width == np.array([120, 120, 100, 9, 103, 101, 120, 120, 104, 11, 102])
    ).all()

    assert (out.strand == np.array([1, -1, 0, 0, 1, -1, 0, 1, -1, -1, 0])).all()


def test_gaps_with_both():
    assert gr is not None

    out = gr.gaps(start=5, end=120)

    assert out is not None
    assert out.seqnames == [
        "chr1",
        "chr1",
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
    ]
    assert (out.start == np.array([5, 5, 5, 112, 5, 5, 5, 5, 5, 110, 5])).all()
    assert (
        out.width == np.array([116, 116, 96, 9, 99, 97, 116, 116, 100, 11, 98])
    ).all()

    assert (out.strand == np.array([1, -1, 0, 0, 1, -1, 0, 1, -1, -1, 0])).all()


def test_disjoin():
    assert gr is not None

    out = gr.disjoin()

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([101, 104, 102, 105, 103])).all()
    assert (out.width == np.array([11, 30, 21, 5, 25])).all()
    assert (out.strand == np.array([0, 1, -1, -1, 0])).all()
