import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

df_gr = pd.DataFrame(
    {
        "seqname": [
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
        "start": range(100, 110),
        "end": range(110, 120),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

gr = GenomicRanges.fromPandas(df_gr)


def test_granges():
    assert gr is not None
    assert gr.len() == df_gr.shape[0]
    assert len(gr) == gr.len()
    assert gr.mcols() is not None
    assert gr.mcols().shape[0] == df_gr.shape[0]
    assert gr.granges() is not None

    test_gr = GenomicRanges.fromPandas(
        pd.DataFrame(
            {
                "seqname": ["chr1", "chr2", "chr3"],
                "start": [100, 115, 119],
                "end": [103, 116, 120],
            }
        )
    )

    hits = gr.nearest(test_gr)
    print(hits)
    assert hits is not None
