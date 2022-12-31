from genomicranges.SeqInfo import SeqInfo
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


seq_obj = {
    "seqnames": ["chr1", "chr2", "chr3",],
    "seqlengths": range(100, 103),
    "isCircular": [random() < 0.5 for _ in range(3)],
    "genome": "hg19",
}


def test_create_SeqInfo():

    seq = SeqInfo(seq_obj)

    assert seq is not None
    assert seq.dims == (3, 3)

    assert seq.isCircular is not None
    assert len(seq.isCircular) == 3

    assert seq.genome is not None
    assert seq.genome == "hg19"

    assert seq.seqlengths is not None
    assert isinstance(seq.seqlengths, dict)
