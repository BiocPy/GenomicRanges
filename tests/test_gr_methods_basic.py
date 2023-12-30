import pytest
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges
from biocframe import BiocFrame
from iranges import IRanges
from random import random
import genomicranges

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
    ranges=IRanges(start=range(100, 110), width=range(110, 120)),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    ),
)


def test_granges():
    assert gr is not None
    assert len(gr) == 10
    assert gr.mcols is not None
    assert len(gr.mcols) == 10
    assert gr.ranges is not None


def test_slices():
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

    subset_gr = gr[5:8]

    assert subset_gr is not None
    assert len(subset_gr) == 3
    assert subset_gr.seqnames == ["chr1", "chr3", "chr3"]


def test_gr_empty_subset():
    gre = GenomicRanges.empty()

    assert gre is not None
    assert isinstance(gre, GenomicRanges)
    assert len(gre) == 0

    subset = gre[0:10]


def test_export():
    df = gr.to_pandas()

    assert df is not None
    assert df.shape[0] == len(gr)
