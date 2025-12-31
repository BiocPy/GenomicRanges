import pytest
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import CompressedGenomicRangesList, GenomicRanges

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
    grl = CompressedGenomicRangesList.from_list(lst=[a, b], names=["a", "b"])

    assert isinstance(grl, CompressedGenomicRangesList)
    assert isinstance(grl["a"], GenomicRanges)
    assert grl["b"] is not None
    assert len(grl) == 2


def test_create_grl_should_fail():
    with pytest.raises(Exception):
        CompressedGenomicRangesList.from_list(lst=[a, 2])


def test_empty_grl():
    grl = CompressedGenomicRangesList.empty(n=100)
    assert isinstance(grl, CompressedGenomicRangesList)

    subset = grl[0:10]
    assert isinstance(subset, CompressedGenomicRangesList)
    assert len(subset) == 10

    subset = grl[1]
    assert isinstance(subset, GenomicRanges)
    assert len(subset) == 0
