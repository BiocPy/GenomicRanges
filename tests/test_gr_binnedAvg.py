import pandas as pd
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

df_src = pd.DataFrame({"seqnames": ["chr1"], "starts": [101], "ends": [109],})

g_src = genomicranges.fromPandas(df_src)

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
        "GC": [random() for _ in range(10)],
    }
)

g_tgt = genomicranges.fromPandas(df_tgt)


def test_binned_average():
    assert g_tgt is not None
    assert g_src is not None

    res = g_tgt.binnedAverage(bins=g_src, scorename="score", outname="binned_score")

    assert res is not None
    assert res.column("binned_score") == [2.6]
