import pytest
from genomicranges import GenomicRanges
from biocframe import BiocFrame
from iranges import IRanges
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


def test_to_pandas_with_names_and_mcols_is_positional():
    # Regression test for a bug where to_pandas() set the index to `names`
    # before concatenating `mcols`, whose index was still the default
    # positional RangeIndex. pandas.concat(axis=1) then outer-joined on the
    # two disjoint indices instead of binding columns positionally, doubling
    # the row count and NaN-splitting every row.
    ranges = IRanges(start=[0, 10, 20], width=[5, 5, 5])
    mcols = BiocFrame({"gene_id": ["g1", "g2", "g3"], "gene_name": ["A", "B", "C"]})
    gr = GenomicRanges(
        seqnames=["1", "1", "1"],
        ranges=ranges,
        strand=["+", "+", "-"],
        names=["g1", "g2", "g3"],
        mcols=mcols,
    )

    df = gr.to_pandas()

    assert df.shape == (3, 7)
    assert list(df.index) == ["g1", "g2", "g3"]
    assert df["gene_id"].notna().all()
    assert df["seqnames"].notna().all()
    assert list(df["gene_id"]) == ["g1", "g2", "g3"]
    assert list(df["gene_name"]) == ["A", "B", "C"]
