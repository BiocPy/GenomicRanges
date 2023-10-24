from genomicranges.Seqinfo import Seqinfo
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


seq_obj = {}


def test_create_Seqinfo():
    circ = [random() < 0.5 for _ in range(3)]
    seq = Seqinfo(
        seqnames=["chr1", "chr2", "chr3"],
        seqlengths=range(100, 103),
        is_circular=circ,
        genome="hg19",
    )

    assert len(seq) == 3
    assert seq.is_circular() == circ
    assert seq.genome() == ["hg19"] * 3
    assert seq.seqlengths() == [100, 101, 102]

    assert seq.seqnames() == ["chr1", "chr2", "chr3"]
    assert seq.seqlengths(as_dict=True) == {"chr1": 100, "chr2": 101, "chr3": 102}
    seq2 = seq.set_seqlengths({"chr2": 500, "chr1": 123, "chr3": 99})
    assert seq2.seqlengths() == [123, 500, 99]

    seq2 = seq.set_is_circular(False)
    assert seq2.is_circular() == [False] * 3
    assert seq2.is_circular(as_dict=True) == {
        "chr1": False,
        "chr2": False,
        "chr3": False,
    }

    seq2 = seq.set_genome(None)
    assert seq2.genome() == [None] * 3
    assert seq2.genome(as_dict=True) == {"chr1": None, "chr2": None, "chr3": None}
