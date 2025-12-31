from random import random

from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import CompressedGenomicRangesList, GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

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


def test_split():
    assert subject is not None

    splits = subject.split(
        [
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
        ]
    )

    assert splits is not None
    assert isinstance(splits, CompressedGenomicRangesList)
    assert len(splits) == 3


def test_to_granges():
    assert subject is not None

    splits = subject.split(
        [
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
        ]
    )

    roundtrip = splits._unlist_data

    assert roundtrip is not None
    assert isinstance(roundtrip, GenomicRanges)
    assert len(roundtrip) == len(subject)
