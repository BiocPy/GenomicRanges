import biocutils as ut
import pytest
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import GenomicRanges, GenomicRangesList

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


def test_create_grl():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    assert isinstance(grl, GenomicRangesList)
    assert isinstance(grl["a"], GenomicRanges)
    assert grl["b"] is not None
    assert len(grl) == 2


def test_create_grl_should_fail():
    with pytest.raises(Exception):
        GenomicRangesList(ranges=[a, 2])


def test_empty_grl():
    grl = GenomicRangesList.empty(n=100)
    assert isinstance(grl, GenomicRangesList)


def test_grl_set_names():
    grl = GenomicRangesList(ranges=[a, b], names=["a", "b"])

    grl_replace_names = grl.set_names(["aa", "bb"])

    assert grl_replace_names is not None
    assert isinstance(grl_replace_names, GenomicRangesList)
    assert list(grl.get_names()) == ["a", "b"]

    assert list(ut.extract_row_names(grl)) == ["a", "b"]
    assert list(ut.extract_row_names(grl_replace_names)) == ["aa", "bb"]
