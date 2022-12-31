import pandas as pd
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

df_subject = pd.DataFrame(
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
        "starts": range(1, 11),
        "ends": [10] * 10,
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

subject = genomicranges.fromPandas(df_subject)

df_query = pd.DataFrame(
    {"seqnames": ["chr2",], "starts": [4], "ends": [6], "strand": ["+"]}
)

query = genomicranges.fromPandas(df_query)


def test_findOverlaps():
    assert subject is not None
    assert query is not None

    res = subject.findOverlaps(query)

    print(res._data)

    assert res is not None
    assert res.shape[1] - query.shape[1] == 1
    assert res.column("hits") == [[1, 2]]


def test_findOverlaps_queryType():
    assert subject is not None
    assert query is not None

    res = subject.findOverlaps(query, queryType="within")

    print(res._data)

    assert res is not None
    assert res.shape[1] - query.shape[1] == 1
    assert res.column("hits") == [[1, 2]]


def test_countOverlaps():
    assert subject is not None
    assert query is not None

    res = subject.countOverlaps(query)

    assert res is not None
    assert len(res) == 1
    assert res[0] == 2


def test_subsetByOverlaps():

    df_query2 = pd.DataFrame(
        {
            "seqnames": ["chr2", "chr4"],
            "starts": [4, 4],
            "ends": [6, 6],
            "strand": ["+", "-"],
        }
    )

    query2 = genomicranges.fromPandas(df_query2)

    assert subject is not None
    assert query2 is not None

    res = subject.subsetByOverlaps(query2)

    assert res is not None
    print(res.data)
    assert res.shape == (1, 4)
