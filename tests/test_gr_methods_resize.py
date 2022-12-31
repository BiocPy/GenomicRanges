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
        "seqnames": ["chr1", "chr2", "chr3", "chr2", "chr3",],
        "starts": range(101, 106),
        "ends": [112, 123, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

gr = genomicranges.fromPandas(df_gr)


def test_resize_default():
    assert gr is not None

    resized_gr = gr.resize(width=10)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [101, 114, 103, 104, 102]
    assert resized_gr.column("ends") == [110, 123, 112, 113, 111]


def test_resize_fix_end():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end")

    assert resized_gr is not None
    assert resized_gr.column("starts") == [103, 102, 119, 125, 105]
    assert resized_gr.column("ends") == [112, 111, 128, 134, 114]


def test_resize_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, ignoreStrand=True)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [101, 102, 103, 104, 105]
    assert resized_gr.column("ends") == [110, 111, 112, 113, 114]


def test_resize_fix_end_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end", ignoreStrand=True)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [103, 114, 119, 125, 102]
    assert resized_gr.column("ends") == [112, 123, 128, 134, 111]
