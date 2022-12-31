import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

df_src = pd.DataFrame(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3", "chr2"],
        "starts": [101, 102, 103, 104, 109],
        "ends": [112, 103, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

g_src = genomicranges.fromPandas(df_src)


def test_invert():
    assert g_src is not None

    result = g_src.invertStrand()

    assert result is not None
    assert result.shape == g_src.shape


def test_sample():
    assert g_src is not None

    result = g_src.sample(k=3)

    assert result is not None
    assert result.shape[0] == 3


def test_concat():
    df_tgt = pd.DataFrame(
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
            "starts": range(101, 111),
            "ends": range(121, 131),
            "strand": ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
            "score": range(0, 10),
            "GC2": [random() for _ in range(10)],
        }
    )

    g_tgt = genomicranges.fromPandas(df_tgt)
    assert g_tgt is not None

    result = g_src.concat(g_tgt)

    assert result is not None
    assert result.shape[0] == 15
