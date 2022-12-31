import pandas as pd
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


def test_duplicated():
    assert gr is not None

    query_hits = gr.duplicated()

    assert query_hits is not None
    assert query_hits == [
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
    ]


def test_matches():
    assert gr is not None

    query_hits = gr.match(gr[2:5, :])

    assert query_hits is not None
    assert query_hits == [[2], [3], [4]]


def test_isUnsorted():
    assert gr is not None

    result = gr.isUnsorted()

    assert result is True


def test_order():
    assert gr is not None

    order = gr.order()

    assert order == [4, 5, 0, 3, 1, 2, 6, 7, 8, 9]


def test_sort():
    assert gr is not None
    result = gr.sort()

    assert result is not None
    assert result.shape == gr.shape

    result = gr.sort(decreasing=True)

    assert result is not None
    assert result.shape == gr.shape


def test_rank():
    assert gr is not None
    result = gr.rank()

    assert result is not None
    assert result == [2, 4, 5, 3, 0, 1, 6, 7, 8, 9]
