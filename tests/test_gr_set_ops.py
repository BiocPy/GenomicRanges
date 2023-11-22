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

gr0 = GenomicRanges(
    seqnames=["chr1", "chr1"], ranges=IRanges([2, 9], [6, 11]), strand=["+", "-"]
)

gr1 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([5], [6]), strand=["-"])


def test_union():
    gr0 = GenomicRanges(
        seqnames=["chr1", "chr1"], ranges=IRanges([2, 9], [6, 11]), strand=["+", "-"]
    )
    gr1 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([5], [6]), strand=["-"])

    assert gr0 is not None
    assert gr1 is not None

    out = gr0.union(gr1)

    assert out is not None
    assert out.seqnames == ["chr1", "chr1"]
    assert (out.start == np.array([2, 5])).all()
    assert (out.width == np.array([6, 15])).all()
    assert (out.strand == np.array([1, -1])).all()


def test_diff():
    gr0 = GenomicRanges(
        seqnames=["chr1", "chr1"], ranges=IRanges([2, 9], [6, 11]), strand=["+", "-"]
    )

    gr1 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([5], [6]), strand=["-"])

    assert gr0 is not None
    assert gr1 is not None

    out = gr0.setdiff(gr1)

    assert out is not None
    assert out.seqnames == ["chr1", "chr1"]
    assert (out.start == np.array([2, 11])).all()
    assert (out.width == np.array([6, 9])).all()
    assert (out.strand == np.array([1, -1])).all()


def test_intersect():
    gr0 = GenomicRanges(
        seqnames=["chr1", "chr1"], ranges=IRanges([2, 9], [6, 11]), strand=["+", "-"]
    )

    gr1 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([5], [6]), strand=["-"])
    assert gr0 is not None
    assert gr1 is not None

    out = gr0.intersect(gr1)

    assert out is not None
    assert out.seqnames == ["chr1"]
    assert (out.start == np.array([9])).all()
    assert (out.width == np.array([2])).all()
    assert (out.strand == np.array([-1])).all()
