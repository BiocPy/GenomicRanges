import pytest
import pandas as pd
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


def test_nearest():
    assert gr is not None

    test_gr = genomicranges.fromPandas(
        pd.DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [200, 105, 1190],
                "ends": [203, 106, 1200],
            }
        )
    )

    query_hits = gr.nearest(test_gr)

    assert query_hits is not None
    assert query_hits.column("hits") == [[4], [3], []]


def test_precede():
    assert gr is not None

    test_gr = genomicranges.fromPandas(
        pd.DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [200, 105, 10],
                "ends": [203, 106, 20],
            }
        )
    )

    query_hits = gr.precede(test_gr)

    print(query_hits.data)

    assert query_hits is not None
    assert query_hits.column("hits") == [[4], [3], []]


def test_follow():
    assert gr is not None

    test_gr = genomicranges.fromPandas(
        pd.DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [200, 105, 1190],
                "ends": [203, 106, 1200],
            }
        )
    )

    query_hits = gr.follow(test_gr)

    print(query_hits.data)

    assert query_hits is not None
    assert query_hits.column("hits") == [[], [3], []]
