import numpy as np
from iranges import IRanges

from genomicranges import GenomicRanges, CompressedGenomicRangesList

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

x = GenomicRanges(seqnames=["chr1", "chr1", "chrX"], ranges=IRanges([1, 40, 1], [50, 71, 500]), strand=["*", "*", "*"])
y = GenomicRanges(seqnames=["chr1", "chr1"], ranges=IRanges([21, 38], [5, 113]), strand=["*", "*"])


def test_subtract():
    out = x.subtract(y)

    assert out is not None
    assert isinstance(out, CompressedGenomicRangesList)
    assert len(out) == 3

    assert out[0].get_seqnames() == ["chr1", "chr1"]
    assert np.all(out[0]._ranges._start == [1, 26])
    assert np.all(out[0]._ranges.get_end() == [20, 37])
    assert np.all(out[2].get_strand() == [0, 0])

    assert len(out[1]) == 0

    assert out[2].get_seqnames() == ["chrX"]
    assert np.all(out[2]._ranges._start == [1])
    assert np.all(out[2]._ranges.get_end() == [500])
    assert np.all(out[2].get_strand() == [0])

    assert len(out) == len(x)
