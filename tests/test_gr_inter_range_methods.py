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
        "seqnames": ["chr1", "chr2", "chr1", "chr3", "chr2"],
        "starts": [101, 102, 103, 104, 109],
        "ends": [112, 103, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

gr = genomicranges.fromPandas(df_gr)


def test_reduce():
    assert gr is not None

    reduced_gr = gr.reduce()

    assert reduced_gr is not None
    assert reduced_gr.shape == (4, 4)
    assert reduced_gr.column("seqnames") == ["chr1", "chr2", "chr2", "chr3"]
    assert reduced_gr.column("starts") == [101, 102, 109, 104]
    assert reduced_gr.column("ends") == [128, 103, 111, 134]
    assert reduced_gr.column("strand") == ["*", "-", "-", "+"]


def test_reduce_with_gapwidth():
    assert gr is not None

    reduced_gr = gr.reduce(minGapwidth=10)

    assert reduced_gr is not None
    assert reduced_gr.shape == (3, 4)
    assert reduced_gr.column("seqnames") == ["chr1", "chr2", "chr3"]
    assert reduced_gr.column("starts") == [101, 102, 104]
    assert reduced_gr.column("ends") == [128, 111, 134]
    assert reduced_gr.column("strand") == ["*", "-", "+"]


def test_reduce_with_gapwidth_withrevmap():
    assert gr is not None

    reduced_gr = gr.reduce(minGapwidth=10, withRevMap=True)

    assert reduced_gr is not None
    assert reduced_gr.shape == (3, 5)
    assert reduced_gr.column("seqnames") == ["chr1", "chr2", "chr3"]
    assert reduced_gr.column("starts") == [101, 102, 104]
    assert reduced_gr.column("ends") == [128, 111, 134]
    assert reduced_gr.column("strand") == ["*", "-", "+"]
    assert reduced_gr.column("revmap") == [[0, 2], [1, 4], [3]]


def test_range():
    assert gr is not None

    range_gr = gr.range()

    assert range_gr is not None
    assert range_gr.shape == (3, 4)
    assert range_gr.column("seqnames") == ["chr1", "chr2", "chr3"]
    assert range_gr.column("starts") == [101, 102, 104]
    assert range_gr.column("ends") == [128, 111, 134]
    assert range_gr.column("strand") == ["*", "-", "+"]


def test_gaps():
    assert gr is not None

    gapped_gr = gr.gaps()

    assert gapped_gr is not None
    assert gapped_gr.column("seqnames") == ["chr1", "chr2", "chr2", "chr3"]
    assert gapped_gr.column("starts") == [1, 1, 104, 1]
    assert gapped_gr.column("ends") == [100, 101, 108, 103]
    assert gapped_gr.column("strand") == ["*", "-", "-", "+"]


def test_gaps_with_start():
    assert gr is not None

    gapped_gr = gr.gaps(start=5)

    assert gapped_gr is not None
    assert gapped_gr.column("seqnames") == ["chr1", "chr2", "chr2", "chr3"]
    assert gapped_gr.column("starts") == [5, 5, 104, 5]
    assert gapped_gr.column("ends") == [100, 101, 108, 103]
    assert gapped_gr.column("strand") == ["*", "-", "-", "+"]


def test_gaps_with_start_filter():
    assert gr is not None

    gapped_gr = gr.gaps(start=103)

    assert gapped_gr is not None
    assert gapped_gr.column("seqnames") == ["chr2", "chr3"]
    assert gapped_gr.column("starts") == [104, 103]
    assert gapped_gr.column("ends") == [108, 103]
    assert gapped_gr.column("strand") == ["-", "+"]


def test_gaps_with_end():
    assert gr is not None

    gapped_gr = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

    assert gapped_gr is not None
    assert gapped_gr.column("seqnames") == [
        "chr1",
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr2",
        "chr2",
        "chr2",
        "chr3",
        "chr3",
        "chr3",
    ]
    assert gapped_gr.column("starts") == [1, 1, 1, 1, 1, 1, 104, 112, 1, 1, 1]
    assert gapped_gr.column("ends") == [
        100,
        120,
        120,
        120,
        120,
        101,
        108,
        120,
        120,
        103,
        120,
    ]

    assert gapped_gr.column("strand") == [
        "*",
        "+",
        "-",
        "*",
        "+",
        "-",
        "-",
        "-",
        "*",
        "+",
        "-",
    ]


def test_disjoin():
    assert gr is not None

    disjoin_gr = gr.disjoin()

    assert disjoin_gr is not None
    assert disjoin_gr.column("seqnames") == [
        "chr1",
        "chr1",
        "chr1",
        "chr2",
        "chr2",
        "chr3",
    ]
    assert disjoin_gr.column("starts") == [101, 103, 113, 102, 109, 104]
    assert disjoin_gr.column("ends") == [102, 112, 128, 103, 111, 134]
    assert disjoin_gr.column("strand") == ["*", "*", "*", "-", "-", "+"]
