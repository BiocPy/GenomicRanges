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


def test_shift():
    assert gr is not None

    shifted_gr = gr.shift(shift=10)

    assert shifted_gr is not None
    assert (shifted_gr.start == np.array([111, 112, 113, 114, 115])).all()
    assert (shifted_gr.width == np.array([11, 21, 25, 30, 5])).all()


def test_promoters():
    assert gr is not None

    prom_gr = gr.promoters()

    assert prom_gr is not None
    assert (prom_gr.start == np.array([-1899, -77, -1897, -1896, -90])).all()
    assert (prom_gr.width == np.array([2200] * 5)).all()


def test_restrict():
    assert gr is not None

    restrict_gr = gr.restrict(start=114, end=140)

    assert restrict_gr is not None
    assert (restrict_gr.start == np.array([114] * 3)).all()
    assert (restrict_gr.width == np.array([9, 14, 20])).all()

    restrict_gr = gr.restrict(start=114, end=140, keep_all_ranges=True)

    assert restrict_gr is not None
    assert (restrict_gr.start == np.array([114] * 5)).all()
    assert (restrict_gr.width == np.array([0, 9, 14, 20, 0])).all()

    restrict_gr = gr.restrict(start=1200)

    assert restrict_gr is not None
    assert len(restrict_gr) == 0


def test_narrow():
    assert gr is not None

    narrow_gr = gr.narrow(start=2, end=3)

    assert narrow_gr is not None
    assert (narrow_gr.start == np.array([102, 103, 104, 105, 106])).all()
    assert (narrow_gr.width == np.array([2] * 5)).all()

    narrow_gr = gr.narrow(start=2)

    assert narrow_gr is not None
    assert (narrow_gr.start == np.array([102, 103, 104, 105, 106])).all()
    assert (narrow_gr.width == np.array([10, 20, 24, 29, 4])).all()

    narrow_gr = gr.narrow(start=2, width=3)

    assert narrow_gr is not None
    assert (narrow_gr.start == np.array([102, 103, 104, 105, 106])).all()
    assert (narrow_gr.width == np.array([3] * 5)).all()

    narrow_gr = gr.narrow(end=2)

    assert narrow_gr is not None
    assert (narrow_gr.start == np.array([101, 102, 103, 104, 105])).all()
    assert (narrow_gr.width == np.array([2] * 5)).all()

    narrow_gr = gr.narrow(end=4, width=3)

    assert narrow_gr is not None
    assert (narrow_gr.start == np.array([102, 103, 104, 105, 106])).all()
    assert (narrow_gr.width == np.array([3] * 5)).all()
