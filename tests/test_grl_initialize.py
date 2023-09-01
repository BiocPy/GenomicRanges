import pytest
import pandas as pd
from genomicranges import GenomicRanges, GenomicRangesList
from biocframe import BiocFrame
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

a = GenomicRanges(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3"],
        "starts": [1, 3, 2, 4],
        "ends": [10, 30, 50, 60],
        "strand": ["-", "+", "*", "+"],
        "score": [1, 2, 3, 4],
    }
)

b = GenomicRanges(
    {
        "seqnames": ["chr2", "chr4", "chr5"],
        "starts": [3, 6, 4],
        "ends": [30, 50, 60],
        "strand": ["-", "+", "*"],
        "score": [2, 3, 4],
    }
)


def test_create_grl():
    grl = GenomicRangesList(a=a, b=b)

    assert isinstance(grl, GenomicRangesList)
    assert isinstance(grl["a"], GenomicRanges)
    assert grl["b"] is not None
    assert len(grl) == 2


def test_create_grl_should_fail():
    with pytest.raises(Exception):
        GenomicRangesList(a=a, b=2)
