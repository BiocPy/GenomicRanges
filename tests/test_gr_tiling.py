import pytest
import pandas as pd
from genomicranges import GenomicRanges
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


def test_tile():
    assert gr is not None

    tiles = gr.tile(n=2)
    assert tiles is not None
    assert isinstance(tiles, GenomicRanges)
    assert len(tiles) == 10


def test_slide():
    assert gr is not None

    tiles = gr.sliding_windows(width=3)

    assert tiles is not None
    assert isinstance(tiles, GenomicRanges)
    assert len(tiles) == 82


def test_tile_genome():
    seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

    tiles = GenomicRanges.tile_genome(seqlengths=seqlengths, n=10)

    assert tiles is not None
    assert isinstance(tiles, GenomicRanges)
    assert len(tiles) == 30
