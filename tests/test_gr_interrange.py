from random import random

import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges.GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

gr = GenomicRanges(
    seqnames=["chr1", "chr2", "chr2", "chr2", "chr1", "chr1", "chr3", "chr3", "chr3", "chr3"],
    ranges=IRanges([x for x in range(1, 11)], [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(1, 11),
            "GC": [random() for _ in range(10)],
        }
    ),
)


def test_reduce():
    assert gr is not None

    out = gr.reduce()

    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([6, 1, 5, 2, 4, 7, 9])).all()
    assert (out.end == np.array([10] * 7)).all()
    assert (out.strand == np.array([1, -1, 0, 1, 0, 1, -1])).all()

    out = gr.reduce(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr3"]
    assert (out.start == np.array([1, 2, 7])).all()
    assert (out.end == np.array([10] * 3)).all()
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
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([6, 1, 5, 2, 4, 7, 9])).all()
    assert (out.end == np.array([10] * 7)).all()
    assert (out.strand == np.array([1, -1, 0, 1, 0, 1, -1])).all()


def test_reduce_with_gapwidth_with_reverse_map():
    assert gr is not None

    out = gr.reduce(min_gap_width=10, with_reverse_map=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([6, 1, 5, 2, 4, 7, 9])).all()
    assert (out.end == np.array([10] * 7)).all()
    assert (out.strand == np.array([1, -1, 0, 1, 0, 1, -1])).all()
    assert out.mcols.column("revmap") == [[5], [0], [4], [1, 2], [3], [6, 7], [8, 9]]


def test_range():
    assert gr is not None

    out = gr.range()

    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([6, 1, 5, 2, 4, 7, 9])).all()
    assert (out.end == np.array([10] * 7)).all()
    assert (out.strand == np.array([1, -1, 0, 1, 0, 1, -1])).all()

    out = gr.range(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr2", "chr3"]
    assert (out.start == np.array([1, 2, 7])).all()
    assert (out.end == np.array([10] * 3)).all()
    assert (out.strand == np.array([0, 0, 0])).all()


def test_gaps():
    assert gr is not None

    out = gr.gaps()
    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert (out.start == np.array([1] * 6)).all()
    assert (out.end == np.array([5, 4, 1, 3, 6, 8])).all()
    assert (out.strand == np.array([1, 0, 1, 0, 1, -1])).all()

    out = gr.gaps(ignore_strand=True)
    assert out is not None
    assert out.seqnames == ["chr2", "chr3"]
    assert (out.start == np.array([1, 1])).all()
    assert (out.end == np.array([1, 6])).all()
    assert (out.strand == np.array([0, 0])).all()


def test_gaps_with_start():
    assert gr is not None

    out = gr.gaps(start=5)

    assert out is not None
    assert out.seqnames == ["chr1", "chr3", "chr3"]
    assert (out.start == np.array([5] * 3)).all()
    assert (out.end == np.array([5, 6, 8])).all()
    assert (out.strand == np.array([1, 1, -1])).all()


def test_gaps_with_start_filter():
    assert gr is not None

    out = gr.gaps(start=103)

    assert out is not None
    assert len(out) == 0


def test_gaps_with_end():
    assert gr is not None

    out = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

    assert out is not None
    assert out.seqnames == [
        "chr1",
        "chr1",
        "chr1",
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr2",
        "chr2",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
        "chr3",
    ]
    assert (out.start == np.array([1, 11, 11, 1, 11, 1, 11, 1, 1, 11, 1, 11, 1, 11, 1])).all()
    assert (out.width == np.array([5, 110, 110, 4, 110, 1, 110, 120, 3, 110, 6, 110, 8, 110, 120])).all()
    assert (out.strand == np.array([1, 1, -1, 0, 0, 1, 1, -1, 0, 0, 1, 1, -1, -1, 0])).all()


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
        "chr3",
    ]
    assert (out.start == np.array([5, 11, 11, 11, 11, 5, 11, 5, 11, 5, 11, 5])).all()
    assert (out.width == np.array([1, 110, 110, 110, 110, 116, 110, 2, 110, 4, 110, 116])).all()
    assert (out.strand == np.array([1, 1, -1, 0, 1, -1, 0, 1, 1, -1, -1, 0])).all()


def test_disjoin():
    assert gr is not None

    out = gr.disjoin()

    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3", "chr3", "chr3", "chr3"]
    assert (out.start == np.array([6, 1, 5, 2, 3, 4, 7, 8, 9, 10])).all()
    assert (out.end == np.array([10, 10, 10, 2, 10, 10, 7, 10, 9, 10])).all()
    assert (out.strand == np.array([1, -1, 0, 1, 1, 0, 1, 1, -1, -1])).all()

    out = gr.disjoin(ignore_strand=True)

    assert out is not None
    assert out.seqnames == ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3", "chr3", "chr3", "chr3"]
    assert (out.start == np.array([1, 5, 6, 2, 3, 4, 7, 8, 9, 10])).all()
    assert (out.end == np.array([4, 5, 10, 2, 3, 10, 7, 8, 9, 10])).all()
    assert (out.strand == np.array([0] * 10)).all()


def test_is_disjoint():
    assert gr is not None

    out = gr.is_disjoint()
    assert out is False

    out = gr.is_disjoint(ignore_strand=True)
    assert out is False


def test_disjoint_bins():
    assert gr is not None

    out = gr.disjoint_bins()
    assert np.all(out == [0, 0, 1, 0, 0, 0, 0, 1, 0, 1])

    out = gr.disjoint_bins(ignore_strand=True)
    assert np.all(out == [0, 0, 1, 2, 1, 2, 0, 1, 2, 3])
