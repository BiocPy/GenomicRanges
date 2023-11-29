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

query = GenomicRanges(seqnames=["chr1"], ranges=IRanges([101], [109]))

subject = GenomicRanges(
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
    ranges=IRanges(range(101, 111), range(121, 131)),
    strand=["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    ),
)


def test_binned_average():
    assert query is not None
    assert subject is not None

    res = subject.binned_average(bins=query, scorename="score", outname="binned_score")

    assert res is not None
    assert res.mcols.get_column("binned_score") == [2]
