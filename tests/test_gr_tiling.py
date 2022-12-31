import pandas as pd
import genomicranges
from random import random

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


def test_tile():
    assert g_src is not None

    tiles = g_src.tile(n=2)

    assert tiles is not None
    assert tiles.shape == (10, 4)


def test_slide():
    assert g_src is not None

    tiles = g_src.slidingWindows(width=10)

    assert tiles is not None
    assert tiles.shape == (44, 4)


def test_tileGenome():
    seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

    tiles = genomicranges.tileGenome(seqlengths=seqlengths, n=10)

    assert tiles is not None
    assert tiles.shape == (30, 4)
