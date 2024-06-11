import pytest
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random
import polars as pl

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_from_polars():
    df_src = pl.DataFrame(
        {
            "seqnames": ["chr1"],
            "starts": [101],
            "widths": [8],
        }
    )

    g_src = GenomicRanges.from_polars(df_src)

    assert g_src is not None
    assert isinstance(g_src, GenomicRanges)
    assert g_src.mcols is not None
    assert isinstance(g_src.mcols, BiocFrame)
    assert len(g_src) == 1
    assert g_src.names is None
    assert g_src.strand is not None


def test_from_polars_should_fail():
    df_gr = pl.DataFrame(
        {
            "starts": range(100, 110),
            "ends": range(110, 120),
            "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    )
    with pytest.raises(Exception):
        GenomicRanges.from_polars(df_gr)


def test_to_polars():
    df_src = pl.DataFrame(
        {
            "seqnames": ["chr1"],
            "starts": [101],
            "widths": [8],
        }
    )

    g_src = GenomicRanges.from_polars(df_src)

    roundtrip = g_src.to_polars()

    assert len(set(roundtrip.columns).difference(df_src.columns)) == 2
    assert len(roundtrip) == len(df_src)


def test_to_polars_complex():
    gr = GenomicRanges(
        seqnames=[
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
        ranges=IRanges(start=range(100, 110), width=range(110, 120)),
        strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        mcols=BiocFrame(
            {
                "score": range(0, 10),
                "GC": [random() for _ in range(10)],
            }
        ),
    )

    roundtrip = gr.to_polars()

    assert len(roundtrip.columns) == 7
    assert len(roundtrip) == len(gr)
