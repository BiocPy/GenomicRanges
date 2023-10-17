import pytest
import pandas as pd
from genomicranges import GenomicRanges, GenomicRangesList
from biocgenerics.combine import combine
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


def test_is_empty_False():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False


def test_is_empty_slice():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False

    sgrl = grl[0:1, :]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRangesList)
    assert len(sgrl) == 1


def test_is_empty_True():
    grl = GenomicRangesList(GenomicRanges.empty(), range_lengths=[0] * 10)

    assert grl.is_empty() == True
    assert len(grl) == 10


def test_is_empty_True_slice():
    grl = GenomicRangesList(GenomicRanges.empty(), range_lengths=[0] * 10)

    sgrl = grl[1:5]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRangesList)
    assert len(sgrl) == 4


def test_nrows():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    nrows = grl.element_nrows()
    assert isinstance(nrows, dict)
    assert list(nrows.keys()) == ["a", "b"]
    assert list(nrows.values()) == [4, 3]


def test_props():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    props = ["start", "end", "strand", "genome", "score"]

    for prop in props:
        v = getattr(grl, prop)

        assert isinstance(v, dict)

    assert isinstance(grl.mcols, BiocFrame)


def test_to_pandas():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    df = grl.to_pandas()
    assert len(df) == 7
    assert df.index.tolist() == [
        "a",
        "a",
        "a",
        "a",
        "b",
        "b",
        "b",
    ]


def test_combine():
    grla = GenomicRangesList(ranges=[a], names=["a"])
    grlb = GenomicRangesList(ranges=[b, a], names=["b", "c"])

    grlc = grla.combine(grlb)

    cgrl = combine(grla, grlb)

    assert len(grlc) == len(cgrl) == 3
