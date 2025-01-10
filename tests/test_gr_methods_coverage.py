from genomicranges.GenomicRanges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random
import numpy as np

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


def test_coverage_default():
    assert gr is not None

    res = gr.coverage()

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 10
    assert np.all(res["chr1"] == [1, 1, 1, 1, 2, 3, 3, 3, 3, 3])
    assert len(res["chr2"]) == 10
    assert np.all(res["chr2"] == [0, 1, 2, 3, 3, 3, 3, 3, 3, 3])
    assert len(res["chr3"]) == 10
    assert np.all(res["chr3"] == [0, 0, 0, 0, 0, 0, 1, 2, 3, 4])


def test_coverage_shift():
    assert gr is not None

    res = gr.coverage(shift=10)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 20
    assert len(res["chr2"]) == 20
    assert len(res["chr3"]) == 20


def test_coverage_shift_width():
    assert gr is not None

    res = gr.coverage(shift=10, width=5)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 5
    assert len(res["chr2"]) == 5
    assert len(res["chr3"]) == 5
