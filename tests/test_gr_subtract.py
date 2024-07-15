import pytest
from genomicranges import GenomicRanges, GenomicRangesList
from iranges import IRanges
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

x = GenomicRanges(
    seqnames=["chr1", "chr1"], ranges=IRanges([2, 9], [6, 11]), strand=["+", "-"]
)
y = GenomicRanges(seqnames=["chr1"], ranges=IRanges([5], [6]), strand=["-"])


def test_subtract():
    out = x.subtract(y)

    assert out is not None
    assert isinstance(out, GenomicRangesList)
    assert (out[0].start == np.array([2])).all()
    assert (out[1].start == np.array([11])).all()

    assert (out[0].width == np.array([6])).all()
    assert (out[1].width == np.array([9])).all()

    assert len(out) == len(x)
