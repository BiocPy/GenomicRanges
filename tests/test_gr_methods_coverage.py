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


def test_coverage_default():
    assert gr is not None

    res = gr.coverage()

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 128
    assert len(res["chr2"]) == 111
    assert len(res["chr3"]) == 134


def test_coverage_shift():
    assert gr is not None

    res = gr.coverage(shift=10)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 138
    assert len(res["chr2"]) == 121
    assert len(res["chr3"]) == 144


def test_coverage_shift_width():
    assert gr is not None

    res = gr.coverage(shift=10, width=5)

    assert res is not None
    assert len(res.keys()) == 3
    assert len(res["chr1"]) == 5
    assert len(res["chr2"]) == 5
    assert len(res["chr3"]) == 5
