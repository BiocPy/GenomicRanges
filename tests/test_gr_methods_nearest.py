import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

df_gr = pd.DataFrame(
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

gr = GenomicRanges.fromPandas(df_gr)


def test_granges():
    assert gr is not None

    test_gr = GenomicRanges.fromPandas(
        pd.DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [100, 115, 119],
                "ends": [103, 116, 120],
            }
        )
    )

    hits = gr.nearest(test_gr)
    assert hits is not None
