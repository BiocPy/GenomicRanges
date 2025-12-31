import biocutils as ut
import numpy as np
import pytest
from biocframe import BiocFrame
from compressed_lists import CompressedCharacterList, CompressedList, CompressedNumpyList
from iranges import CompressedIRangesList, IRanges

from genomicranges import CompressedGenomicRangesList, GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


@pytest.fixture
def grl_named():
    gr1 = GenomicRanges(
        seqnames=["chr1", "chr2"],
        ranges=IRanges([1, 10], [5, 5]),
        strand=["+", "-"],
        mcols=BiocFrame({"score": [1, 2]}),
    )
    gr2 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([20], [5]), strand=["*"], mcols=BiocFrame({"score": [3]}))
    return CompressedGenomicRangesList.from_list([gr1, gr2], names=["a", "b"])


@pytest.fixture
def grl_unnamed():
    gr1 = GenomicRanges(seqnames=["chr1"], ranges=IRanges([1], [10]), strand=["+"])
    gr2 = GenomicRanges(seqnames=["chr2"], ranges=IRanges([5], [10]), strand=["-"])
    return CompressedGenomicRangesList.from_list([gr1, gr2])


def test_grl_seqnames(grl_named):
    seqs = grl_named.seqnames
    assert isinstance(seqs, CompressedCharacterList)
    assert len(seqs) == 2
    assert seqs[0] == ut.StringList(["chr1", "chr2"])
    assert seqs[1] == ut.StringList(["chr1"])

    assert isinstance(seqs.unlist(), ut.StringList)
    assert seqs.unlist()[0] == "chr1"


def test_grl_strand(grl_named):
    strands = grl_named.strand
    assert isinstance(strands, CompressedCharacterList)
    assert len(strands) == 2
    assert strands[0] == ut.StringList(["+", "-"])
    assert strands[1] == ut.StringList(["*"])


def test_grl_ranges(grl_named):
    rngs = grl_named.ranges
    assert isinstance(rngs, CompressedIRangesList)
    assert len(rngs) == 2

    assert isinstance(rngs[0], IRanges)
    assert np.allclose(rngs[0].start, np.array([1, 10]))
    assert np.allclose(rngs[0].width, np.array([5, 5]))


def test_grl_coordinates(grl_named):
    starts = grl_named.start
    assert isinstance(starts, CompressedNumpyList)
    assert np.allclose(starts[0], np.asarray([1, 10]))
    assert np.allclose(starts[1], np.asarray([20]))

    ends = grl_named.end
    assert isinstance(ends, CompressedNumpyList)
    assert np.allclose(ends[0], np.asarray([5, 14]))
    assert np.allclose(ends[1], np.asarray([24]))

    widths = grl_named.width
    assert isinstance(widths, CompressedList)
    assert np.allclose(widths[0], np.asarray([5, 5]))
    assert np.allclose(widths[1], np.asarray([5]))


def test_grl_seqinfo(grl_named):
    si = grl_named.seqinfo
    assert si is not None
    assert set(si.get_seqnames()) == {"chr1", "chr2"}


def test_grl_stack_named(grl_named):
    stacked = grl_named.stack(index_column_name="sample_id")

    assert isinstance(stacked, GenomicRanges)
    assert len(stacked) == 3

    assert "sample_id" in stacked.mcols.column_names
    indices = stacked.mcols.get_column("sample_id")
    assert indices.tolist() == ["a", "a", "b"]

    assert "score" in stacked.mcols.column_names
    assert stacked.mcols.get_column("score") == [1, 2, 3]


def test_grl_stack_unnamed(grl_unnamed):
    stacked = grl_unnamed.stack(index_column_name="idx")

    assert isinstance(stacked, GenomicRanges)
    assert len(stacked) == 2

    assert "idx" in stacked.mcols.column_names
    indices = stacked.mcols.get_column("idx")
    assert indices.tolist() == [0, 1]
