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


def test_flank_default():
    assert gr is not None

    flanked_gr = gr.flank(width=10)

    assert flanked_gr is not None
    assert flanked_gr.column("starts") == [91, 124, 93, 94, 112]
    assert flanked_gr.column("ends") == [100, 133, 102, 103, 121]


def test_flank_start_false():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False)

    assert flanked_gr is not None
    assert flanked_gr.column("starts") == [113, 92, 129, 135, 95]
    assert flanked_gr.column("ends") == [122, 101, 138, 144, 104]


def test_flank_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, both=True)

    assert flanked_gr is not None
    assert flanked_gr.column("starts") == [91, 114, 93, 94, 102]
    assert flanked_gr.column("ends") == [110, 133, 112, 113, 121]


def test_flank_start_false_and_both_true():
    assert gr is not None

    flanked_gr = gr.flank(width=10, start=False, both=True)

    assert flanked_gr is not None
    assert flanked_gr.column("starts") == [103, 92, 119, 125, 95]
    assert flanked_gr.column("ends") == [122, 111, 138, 144, 114]
