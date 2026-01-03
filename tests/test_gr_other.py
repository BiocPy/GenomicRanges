from random import random

import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges.GenomicRanges import GenomicRanges
import pytest

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


def test_shift():
    assert gr is not None

    shifted_gr = gr.shift(shift=10)

    assert shifted_gr is not None
    assert isinstance(shifted_gr, GenomicRanges)
    assert shifted_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(shifted_gr.start == np.array([11, 12, 13, 14, 15, 16, 17, 18, 19, 20]))
    assert np.all(shifted_gr.end == np.array([20] * 10))


def test_promoters():
    assert gr is not None

    prom_gr = gr.promoters()

    assert prom_gr is not None
    assert prom_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(prom_gr.start == np.array([-189, -1998, -1997, -1996, -1995, -1994, -1993, -1992, -189, -189]))
    assert np.all(prom_gr.end == np.array([2010, 201, 202, 203, 204, 205, 206, 207, 2010, 2010]))


def test_terminators():
    assert gr is not None

    prom_gr = gr.terminators()

    assert prom_gr is not None
    assert prom_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(prom_gr.start == np.array([-198, -1990, -1990, -1990, -1990, -1990, -1990, -1990, -190, -189]))
    assert np.all(prom_gr.end == np.array([2001, 209, 209, 209, 209, 209, 209, 209, 2009, 2010]))


def test_restrict():
    assert gr is not None

    restrict_gr = gr.restrict(start=3, end=7)

    assert restrict_gr is not None
    assert np.all(restrict_gr.start == [3, 3, 3, 4, 5, 6, 7, 8])
    assert np.all(restrict_gr.end == [7, 7, 7, 7, 7, 7, 7, 7])

    restrict_gr = gr.restrict(start=4, end=8, keep_all_ranges=True)

    assert restrict_gr is not None
    assert np.all(restrict_gr.start == np.array([4, 4, 4, 4, 5, 6, 7, 8, 9, 9]))
    assert np.all(restrict_gr.end == np.array([8] * 10))

    restrict_gr = gr.restrict(start=1200)

    assert restrict_gr is not None
    assert len(restrict_gr) == 0

    starts = {
        "chr1": 4,
        "chr2": 5,
    }

    ends = {"chr2": 8, "chr3": 9}
    restrict_gr = gr.restrict(start=starts, end=ends)
    assert restrict_gr is not None
    assert np.all(restrict_gr.start == np.array([4, 5, 5, 5, 5, 6, 7, 8, 9, 10]))
    assert np.all(restrict_gr.end == np.array([10, 8, 8, 8, 10, 10, 9, 9, 9, 9]))


def test_narrow():
    assert gr is not None

    with pytest.raises(Exception):
        narrow_gr = gr.narrow(start=2, end=3)

    narrow_gr = gr.narrow(start=2)

    assert narrow_gr is not None
    assert np.all(narrow_gr.start == np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11]))
    assert np.all(narrow_gr.end == np.array([10] * 10))

    with pytest.raises(Exception):
        narrow_gr = gr.narrow(start=2, width=3)

    with pytest.raises(Exception):
        narrow_gr = gr.narrow(end=2)

    with pytest.raises(Exception):
        narrow_gr = gr.narrow(end=4, width=3)
