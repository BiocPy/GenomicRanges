import pandas as pd
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


def test_union():
    assert g_src is not None
    assert g_tgt is not None

    union_gr = g_src.union(g_tgt)

    assert union_gr is not None
    assert union_gr.shape == (6, 4)
    assert union_gr.column("seqnames") == [
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr3",
        "chr3",
    ]
    assert union_gr.column("starts") == [101, 106, 104, 102, 104, 109]
    assert union_gr.column("ends") == [128, 126, 124, 123, 134, 130]
    assert union_gr.column("strand") == ["*", "+", "*", "-", "+", "-"]


def test_intersect():
    assert g_src is not None
    assert g_tgt is not None

    int_gr = g_src.intersect(g_tgt)

    assert int_gr is not None
    assert int_gr.shape == (6, 4)
    assert int_gr.column("seqnames") == ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"]
    assert int_gr.column("starts") == [101, 106, 104, 102, 104, 109]
    assert int_gr.column("ends") == [128, 126, 124, 123, 134, 130]
    assert int_gr.column("strand") == ["*", "+", "*", "-", "+", "-"]


def test_diff():
    assert g_src is not None
    assert g_tgt is not None

    diff_gr = g_src.setdiff(g_tgt)

    assert diff_gr is not None
    assert diff_gr.shape == (3, 4)
    assert diff_gr.column("seqnames") == ["chr1", "chr3", "chr3"]
    assert diff_gr.column("starts") == [126, 104, 129]
    assert diff_gr.column("ends") == [128, 106, 134]
    assert diff_gr.column("strand") == ["*", "+", "+"]
