import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from biocframe import BiocFrame
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_create_gr():
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

    gr = genomicranges.from_pandas(df_gr)

    assert gr is not None
    assert gr.dims[0] == df_gr.shape[0]


def test_nested_bframe():
    obj = {
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
        "score": BiocFrame({"scores": range(0, 10)}),
        "GC": [random() for _ in range(10)],
    }

    gr = GenomicRanges(obj)

    assert gr is not None
    assert gr.dims[1] == len(obj.keys())


def test_should_fail():
    with pytest.raises(Exception):
        df_gr = pd.DataFrame(
            {
                "starts": range(100, 110),
                "ends": range(110, 120),
                "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
                "score": range(0, 10),
                "GC": [random() for _ in range(10)],
            }
        )

        genomicranges.from_pandas(df_gr)


def test_gr_empty():
    gr = GenomicRanges(number_of_rows=100)

    assert gr is not None
    assert isinstance(gr, GenomicRanges)
    assert len(gr) == 100
    assert gr.shape == (100, 0)

    gre = GenomicRanges.empty()

    assert gre is not None
    assert isinstance(gre, GenomicRanges)
    assert len(gre) == 0
    assert gre.shape == (0, 0)
