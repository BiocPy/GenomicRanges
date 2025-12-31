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


def test_is_empty_slice():
    grl = CompressedGenomicRangesList.from_list(lst=[a, b], names=["a", "b"])

    sgrl = grl[0:1]
    assert sgrl is not None
    assert isinstance(sgrl, CompressedGenomicRangesList)
    assert len(sgrl) == 1


def test_slice_by_name():
    grl = CompressedGenomicRangesList.from_list(lst=[a, b], names=["a", "b"])

    sgrl = grl["a"]
    assert sgrl is not None
    assert isinstance(sgrl, GenomicRanges)
    assert len(sgrl) == 4
