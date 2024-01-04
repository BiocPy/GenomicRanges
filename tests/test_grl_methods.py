import pytest
import pandas as pd
from genomicranges import GenomicRanges, GenomicRangesList
from biocutils import combine_sequences
from biocframe import BiocFrame
from iranges import IRanges
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

a = GenomicRanges(
    seqnames=["chr1", "chr2", "chr1", "chr3"],
    ranges=IRanges([1, 3, 2, 4], [10, 30, 50, 60]),
    strand=["-", "+", "*", "+"],
    mcols=BiocFrame({"score": [1, 2, 3, 4]}),
)

b = GenomicRanges(
    seqnames=["chr2", "chr4", "chr5"],
    ranges=IRanges([3, 6, 4], [30, 50, 60]),
    strand=["-", "+", "*"],
    mcols=BiocFrame({"score": [2, 3, 4]}),
)


def test_is_empty_False():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False


def test_is_empty_slice():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False

    sgrl = grl[0:1]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRangesList)
    assert len(sgrl) == 1


def test_slice_by_name():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False

    sgrl = grl[["a"]]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRangesList)
    assert len(sgrl) == 1


def test_slice_by_bool():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert grl.is_empty() == False

    sgrl = grl[[True, False]]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRangesList)
    assert len(sgrl) == 1

    with pytest.raises(Exception):
        grl[[False]]


def test_is_empty_True():
    grl = GenomicRangesList(GenomicRanges.empty(), range_lengths=[0])

    assert grl.is_empty() == True
    assert len(grl) == 1


def test_nrows():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    nrows = grl.element_nrows()
    assert isinstance(nrows, dict)
    assert list(nrows.keys()) == ["a", "b"]
    assert list(nrows.values()) == [4, 3]


def test_props():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    props = ["start", "end", "strand"]

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

    cgrl = combine_sequences(grla, grlb)

    assert len(cgrl) == 3


def test_empty_grl_slice():
    grl = GenomicRangesList.empty(n=100)
    assert isinstance(grl, GenomicRangesList)

    subset = grl[0:10]
    assert isinstance(subset, GenomicRangesList)
    assert len(subset) == 10

    subset = grl[[1, 2, 3]]
    assert isinstance(subset, GenomicRangesList)
    assert len(subset) == 3
