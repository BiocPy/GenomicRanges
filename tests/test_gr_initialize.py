from random import random

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_create_gr():
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
        ranges=IRanges(start=range(100, 110), width=range(110, 120)),
        strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        mcols=BiocFrame(
            {
                "score": range(0, 10),
                "GC": [random() for _ in range(10)],
            }
        ),
        names=["a", "b"] * 5,
    )

    assert gr is not None
    assert len(gr) == 10

    assert list(ut.extract_row_names(gr)) == ["a", "b"] * 5


def test_gr_empty():
    gre = GenomicRanges.empty()

    assert gre is not None
    assert isinstance(gre, GenomicRanges)
    assert len(gre) == 0


def test_create_gr_with_seqnames():
    gr = GenomicRanges(
        seqnames=["chr1"] * 10,
        ranges=IRanges(start=range(100, 110), width=range(110, 120)),
    )

    assert gr is not None
    assert gr._seqnames.dtype == np.uint8

    gr16 = GenomicRanges(
        seqnames=[f"chr{i}" for i in range(500)],
        ranges=IRanges(start=range(0, 500), width=range(10, 510)),
    )

    assert gr16 is not None
    assert gr16._seqnames.dtype == np.uint16

    gr32 = GenomicRanges(
        seqnames=[f"chr{i}" for i in range(2**16 + 1)],
        ranges=IRanges(start=range(0, 2**16 + 1), width=range(10, 2**16 + 11)),
    )

    assert gr32 is not None
    assert gr32._seqnames.dtype == np.uint32
