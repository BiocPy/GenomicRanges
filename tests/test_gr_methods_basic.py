import pytest
import pandas as pd
import numpy as np
from genomicranges.GenomicRanges import GenomicRanges
from random import random
import genomicranges

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

gr = genomicranges.fromPandas(df_gr)


def test_granges():
    assert gr is not None
    assert len(gr) == df_gr.shape[0]
    assert gr.mcols() is not None
    assert len(gr.mcols().keys()) == df_gr.shape[1] - 4
    assert gr.granges() is not None


def test_slices():
    subset_gr = gr[5:8]

    assert subset_gr is not None
    assert len(subset_gr) == 3


def test_export():
    df = gr.toPandas()

    assert df is not None
    assert df.shape[0] == len(gr)
