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


def test_flank_default():
    assert gr is not None

    flanked_gr = gr.flank(width=10)

    print(flanked_gr.__repr__())

    assert flanked_gr is not None
    assert (flanked_gr.start == np.array([91, 123, 93, 94, 110])).all()
    assert (flanked_gr.width == np.array([10, 10, 10, 10, 10])).all()


def test_flank_start_false():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False)

    assert flanked_gr is not None
    assert (flanked_gr.start == np.array([112, 92, 128, 134, 95])).all()
    assert (flanked_gr.width == np.array([10, 10, 10, 10, 10])).all()


def test_flank_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, both=True)

    assert flanked_gr is not None
    assert (flanked_gr.start == np.array([91, 113, 93, 94, 100])).all()
    assert (flanked_gr.width == np.array([20, 20, 20, 20, 20])).all()


def test_flank_start_false_and_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False, both=True)

    assert flanked_gr is not None
    assert (flanked_gr.start == np.array([102, 92, 118, 124, 95])).all()
    assert (flanked_gr.width == np.array([20, 20, 20, 20, 20])).all()
