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
        "starts": [101, 102, 103, 104, 106],
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
    assert resized_gr.column("starts") == [103, 102, 119, 125, 106]
    assert resized_gr.column("ends") == [112, 111, 128, 134, 115]


def test_resize_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, ignoreStrand=True)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [101, 102, 103, 104, 106]
    assert resized_gr.column("ends") == [110, 111, 112, 113, 115]


def test_resize_fix_end_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end", ignoreStrand=True)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [103, 114, 119, 125, 102]
    assert resized_gr.column("ends") == [112, 123, 128, 134, 111]


def test_resize_fix_center_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center", ignoreStrand=True)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [102, 108, 111, 114, 104]
    assert resized_gr.column("ends") == [111, 117, 120, 123, 113]


def test_resize_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center")

    assert resized_gr is not None
    assert resized_gr.column("starts") == [102, 108, 111, 114, 104]
    assert resized_gr.column("ends") == [111, 117, 120, 123, 113]


def test_resize_default_width_odd():
    assert gr is not None

    resized_gr = gr.resize(width=11)

    assert resized_gr is not None
    assert resized_gr.column("starts") == [101, 113, 103, 104, 101]
    assert resized_gr.column("ends") == [111, 123, 113, 114, 111]


def test_resize_default_width_odd_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=11, fix="center")

    assert resized_gr is not None
    assert resized_gr.column("starts") == [101, 107, 110, 114, 103]
    assert resized_gr.column("ends") == [111, 117, 120, 124, 113]
