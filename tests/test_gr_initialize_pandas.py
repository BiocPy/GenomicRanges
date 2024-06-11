import pytest
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random
import pandas as pd

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_from_pandas():
    df_src = pd.DataFrame(
        {
            "seqnames": ["chr1"],
            "starts": [101],
            "widths": [8],
        }
    )

    g_src = GenomicRanges.from_pandas(df_src)

    assert g_src is not None
    assert isinstance(g_src, GenomicRanges)
    assert g_src.mcols is not None
    assert isinstance(g_src.mcols, BiocFrame)
    assert len(g_src) == 1
    assert g_src.names is not None
    assert g_src.strand is not None


def test_from_pandas_should_fail():
    df_gr = pd.DataFrame(
        {
            "starts": range(100, 110),
            "ends": range(110, 120),
            "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    )
    with pytest.raises(Exception):
        GenomicRanges.from_pandas(df_gr)
