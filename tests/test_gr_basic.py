from random import random

import biocutils as ut
import pandas as pd
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

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


def test_granges():
    assert gr is not None
    assert len(gr) == 10
    assert gr.mcols is not None
    assert len(gr.mcols) == 10
    assert gr.ranges is not None


def test_slices():
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

    subset_gr = gr[5:8]

    assert subset_gr is not None
    assert len(subset_gr) == 3
    assert subset_gr.seqnames == ["chr1", "chr3", "chr3"]


def test_gr_empty_subset():
    gre = GenomicRanges.empty()

    assert gre is not None
    assert isinstance(gre, GenomicRanges)
    assert len(gre) == 0

    subset = gre[0:10]
    assert subset is not None
    assert len(subset) == 0


def test_export():
    df = gr.to_pandas()

    assert df is not None
    assert df.shape[0] == len(gr)
    assert df["seqnames"].tolist() == [
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
    ]
    assert df["strand"].tolist() == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]


def test_export_pandas_with_mcols():
    df = GenomicRanges(seqnames=["A"], ranges=IRanges(start=[0], width=[10])).reduce().to_pandas()

    assert df is not None
    assert isinstance(df, pd.DataFrame)


def test_combine():
    g_src = GenomicRanges(
        seqnames=["chr1", "chr2", "chr1", "chr3", "chr2"],
        ranges=IRanges(start=[101, 102, 103, 104, 109], width=[112, 103, 128, 134, 111]),
        strand=["*", "-", "*", "+", "-"],
    )

    g_tgt = GenomicRanges(
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
        ranges=IRanges(start=range(101, 111), width=range(121, 131)),
        strand=["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
        mcols=BiocFrame(
            {
                "score": range(0, 10),
                "GC": [random() for _ in range(10)],
            }
        ),
    )
    assert g_src is not None
    assert g_tgt is not None

    out = ut.combine_sequences(g_src, g_tgt)

    assert out is not None
    assert len(out) == 15
    assert len(out.get_mcols().get_column_names()) == 2


def test_combine_diff():
    a = GenomicRanges(["A"], IRanges([0], [10]))
    b = GenomicRanges(["B"], IRanges([5], [15]))

    assert a is not None
    assert b is not None

    out = ut.combine_sequences(a, b)

    assert out is not None
    assert len(out) == 2
    assert len(out.get_mcols().get_column_names()) == 0
    assert out.get_seqnames() == ["A", "B"]
