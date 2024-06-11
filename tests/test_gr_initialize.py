import pytest
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random
import pandas as pd

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
    )

    assert gr is not None
    assert len(gr) == 10


def test_gr_empty():
    gre = GenomicRanges.empty()

    assert gre is not None
    assert isinstance(gre, GenomicRanges)
    assert len(gre) == 0
