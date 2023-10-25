from genomicranges.Seqinfo import Seqinfo, merge_Seqinfo
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
    assert seq2.genome() == [ None ] * 3
    assert seq2.genome(as_dict=True) == { "chr1": None, "chr2": None, "chr3": None}


def test_merge_Seqinfo():
    seq = Seqinfo(
        seqnames = [ "chr1", "chr2", "chr3" ],
        seqlengths = range(100, 103),
        is_circular = [False, True, False],
        genome = "hg19"
    )

    combined = merge_Seqinfo([seq, seq])
    assert combined.seqnames() == seq.seqnames()
    assert combined.seqlengths() == seq.seqlengths()
    assert combined.is_circular() == seq.is_circular()
    assert combined.genome() == seq.genome()

    seq2 = Seqinfo(
        seqnames = [ "chr3", "chr4", "chr5" ],
        seqlengths = range(100, 103),
        is_circular = [False, True, False],
        genome = "hg38"
    )

    combined = merge_Seqinfo([seq, seq2])
    assert combined.seqnames() == [ "chr1", "chr2", "chr3", "chr4", "chr5" ]
    assert combined.seqlengths() == [ 100, 101, None, 101, 102 ]
    assert combined.is_circular() == [ False, True, False, True, False ]
    assert combined.genome() == [ "hg19", "hg19", None, "hg38", "hg38" ]
