import pandas as pd
from genomicranges import GenomicRanges
from biocutils.rownames import rownames
from biocutils.colnames import colnames
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_names():
    gr = GenomicRanges(
        {
            "seqnames": [
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
            "starts": range(100, 110),
            "ends": range(110, 120),
            "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    )

    assert rownames(gr) == None
    assert colnames(gr) == ["seqnames", "starts", "ends", "strand", "score", "GC"]
