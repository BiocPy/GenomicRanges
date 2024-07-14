import pytest
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges, GenomicRangesList
from biocframe import BiocFrame
from iranges import IRanges
from random import random
import genomicranges

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
    assert isinstance(splits, GenomicRangesList)
    assert len(splits) == 3
    print(splits.element_nrows())
    assert sum(splits.get_range_lengths()) == len(subject)


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

    roundtrip = splits.as_genomic_ranges()

    assert roundtrip is not None
    assert isinstance(roundtrip, GenomicRanges)
    assert len(roundtrip) == len(subject)
