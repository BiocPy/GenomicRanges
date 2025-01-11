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


def test_flank_default():
    assert gr is not None

    flanked_gr = gr.flank(width=10)

    assert flanked_gr is not None
    assert isinstance(flanked_gr, GenomicRanges)
    assert flanked_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(flanked_gr.start == np.array([11, -8, -7, -6, -5, -4, -3, -2, 11, 11]))
    assert np.all(flanked_gr.end == np.array([20, 1, 2, 3, 4, 5, 6, 7, 20, 20]))


def test_flank_start_false():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False)

    assert flanked_gr is not None
    assert flanked_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(flanked_gr.start == np.array([-9, 11, 11, 11, 11, 11, 11, 11, -1, 0]))
    assert np.all(flanked_gr.end == np.array([0, 20, 20, 20, 20, 20, 20, 20, 8, 9]))


def test_flank_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, both=True)

    assert flanked_gr is not None
    assert np.all(flanked_gr.start == np.array([1, -8, -7, -6, -5, -4, -3, -2, 1, 1]))
    assert np.all(flanked_gr.end == np.array([20, 11, 12, 13, 14, 15, 16, 17, 20, 20]))


def test_flank_start_false_and_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False, both=True)

    assert flanked_gr is not None
    assert np.all(flanked_gr.start == np.array([-9, 1, 1, 1, 1, 1, 1, 1, -1, 0]))
    assert np.all(flanked_gr.end == np.array([10, 20, 20, 20, 20, 20, 20, 20, 18, 19]))
