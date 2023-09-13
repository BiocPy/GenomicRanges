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

subject = genomicranges.from_pandas(df_subject)

df_query = pd.DataFrame(
    {
        "seqnames": [
            "chr2",
        ],
        "starts": [4],
        "ends": [6],
        "strand": ["+"],
    }
)

query = genomicranges.from_pandas(df_query)


def test_find_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query)

    print(res._data)

    assert res is not None
    assert res.shape[1] - query.shape[1] == 1
    assert res.column("hits") == [[1, 2]]


def test_find_overlaps_query_type():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query, query_type="within")

    print(res._data)

    assert res is not None
    assert res.shape[1] - query.shape[1] == 1
    assert res.column("hits") == [[1, 2]]


def test_count_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.count_overlaps(query)

    assert res is not None
    assert len(res) == 1
    assert res[0] == 2


def test_subset_by_overlaps():
    df_query2 = pd.DataFrame(
        {
            "seqnames": ["chr2", "chr4"],
            "starts": [4, 4],
            "ends": [6, 6],
            "strand": ["+", "-"],
        }
    )

    query2 = genomicranges.from_pandas(df_query2)

    assert subject is not None
    assert query2 is not None

    res = subject.subset_by_overlaps(query2)

    assert res is not None
    assert res.shape == (1, 4)
