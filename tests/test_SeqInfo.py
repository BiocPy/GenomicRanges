from genomicranges.SeqInfo import SeqInfo, merge_SeqInfo
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


seq_obj = {}


def test_create_SeqInfo():
    circ = [random() < 0.5 for _ in range(3)]
    seq = SeqInfo(
        seqnames=["chr1", "chr2", "chr3"],
        seqlengths=range(100, 103),
        is_circular=circ,
        genome="hg19",
    )

    assert len(seq) == 3
    assert seq.get_is_circular() == circ
    assert seq.get_genome() == ["hg19"] * 3
    assert seq.get_seqlengths() == [100, 101, 102]

    assert seq.get_seqnames() == ["chr1", "chr2", "chr3"]
    seq2 = seq.set_seqlengths({"chr2": 500, "chr1": 123, "chr3": 99})
    assert seq2.get_seqlengths() == [123, 500, 99]

    seq2 = seq.set_is_circular(False)
    assert seq2.get_is_circular() == [False] * 3

    seq2 = seq.set_genome(None)
    assert seq2.get_genome() == [None] * 3


def test_merge_SeqInfo():
    seq = SeqInfo(
        seqnames=["chr1", "chr2", "chr3"],
        seqlengths=range(100, 103),
        is_circular=[False, True, False],
        genome="hg19",
    )

    combined = merge_SeqInfo([seq, seq])
    assert combined.get_seqnames() == seq.get_seqnames()
    assert combined.get_seqlengths() == seq.get_seqlengths()
    assert combined.get_is_circular() == seq.get_is_circular()
    assert combined.get_genome() == seq.get_genome()

    seq2 = SeqInfo(
        seqnames=["chr3", "chr4", "chr5"],
        seqlengths=range(100, 103),
        is_circular=[False, True, False],
        genome="hg38",
    )

    combined = merge_SeqInfo([seq, seq2])
    assert combined.get_seqnames() == ["chr1", "chr2", "chr3", "chr4", "chr5"]
    assert combined.get_seqlengths() == [100, 101, None, 101, 102]
    assert combined.get_is_circular() == [False, True, False, True, False]
    assert combined.get_genome() == ["hg19", "hg19", None, "hg38", "hg38"]
