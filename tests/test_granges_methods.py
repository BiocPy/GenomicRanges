import pytest
import pandas as pd
import numpy as np
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


def test_numpy_ufunc():
    res  = np.sin(gr)
    assert res is not None
    assert res.len() == df_gr.shape[0]
    assert len(res) == res.len()
    assert res.mcols() is not None
    assert res.mcols().shape[0] == df_gr.shape[0]
    assert res.granges() is not None

