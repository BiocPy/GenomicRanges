import pytest
from genomicranges.GenomicRanges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random

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
    ranges=IRanges(start=range(100, 110), width=[10] * 10),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    ),
)


def test_coverage_default():
    assert gr is not None

    res = gr.coverage()

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 115
    assert len(res["chr2"]) == 113
    assert len(res["chr3"]) == 119


def test_coverage_shift():
    assert gr is not None

    res = gr.coverage(shift=10)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 125
    assert len(res["chr2"]) == 123
    assert len(res["chr3"]) == 129


def test_coverage_shift_width():
    assert gr is not None

    res = gr.coverage(shift=10, width=5)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 5
    assert len(res["chr2"]) == 5
    assert len(res["chr3"]) == 5
