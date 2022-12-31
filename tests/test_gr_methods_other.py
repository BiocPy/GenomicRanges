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


def test_shift():
    assert gr is not None

    shifted_gr = gr.shift(shift=10)

    assert shifted_gr is not None
    assert shifted_gr.column("starts") == [111, 112, 113, 114, 115]
    assert shifted_gr.column("ends") == [122, 133, 138, 144, 121]


def test_promoters():
    assert gr is not None

    prom_gr = gr.promoters()

    assert prom_gr is not None
    assert prom_gr.column("starts") == [-1899, -76, -1897, -1896, -88]
    assert prom_gr.column("ends") == [300, 2123, 302, 303, 2111]


def test_restrict():
    assert gr is not None

    restrict_gr = gr.restrict(start=114, end=140)

    assert restrict_gr is not None
    assert restrict_gr.column("starts") == [114] * 3
    assert restrict_gr.column("ends") == [123, 128, 134]

    restrict_gr = gr.restrict(start=114, end=140, keepAllRanges=True)

    assert restrict_gr is not None
    assert restrict_gr.column("starts") == [114] * 5
    assert restrict_gr.column("ends") == [112, 123, 128, 134, 111]

    restrict_gr = gr.restrict(start=1200)

    assert restrict_gr is not None
    assert restrict_gr.dims == (0, 6)


def test_narrow():
    assert gr is not None

    narrow_gr = gr.narrow(start=2, end=3)

    assert narrow_gr is not None
    assert narrow_gr.column("starts") == [102, 103, 104, 105, 106]
    assert narrow_gr.column("ends") == [103, 104, 105, 106, 107]

    narrow_gr = gr.narrow(start=2)

    assert narrow_gr is not None
    assert narrow_gr.column("starts") == [102, 103, 104, 105, 106]
    assert narrow_gr.column("ends") == [112, 123, 128, 134, 111]

    narrow_gr = gr.narrow(start=2, width=3)

    assert narrow_gr is not None
    assert narrow_gr.column("starts") == [102, 103, 104, 105, 106]
    assert narrow_gr.column("ends") == [104, 105, 106, 107, 108]

    narrow_gr = gr.narrow(end=2)

    assert narrow_gr is not None
    assert narrow_gr.column("starts") == [101, 102, 103, 104, 105]
    assert narrow_gr.column("ends") == [102, 103, 104, 105, 106]

    narrow_gr = gr.narrow(end=4, width=3)

    assert narrow_gr is not None
    assert narrow_gr.column("starts") == [102, 103, 104, 105, 106]
    assert narrow_gr.column("ends") == [104, 105, 106, 107, 108]
