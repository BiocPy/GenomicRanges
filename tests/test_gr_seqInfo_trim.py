import pandas as pd
from genomicranges.SeqInfo import SeqInfo
from random import random
from genomicranges.GenomicRanges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


gr = GenomicRanges(
    seqnames=[
        "chr1",
        "chr2",
        "chr3",
        "chr2",
        "chr3",
    ],
    ranges=IRanges([x for x in range(101, 106)], [11, 21, 25, 30, 5]),
    strand=["*", "-", "*", "+", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 5),
            "GC": [random() for _ in range(5)],
        }
    ),
)

seq_obj = SeqInfo(
    seqnames=["chr1", "chr2", "chr3"],
    seqlengths=[110, 112, 118],
    is_circular=[True, True, False],
    genome="hg19",
)


def test_gr_seqInfo():
    assert gr is not None
    assert gr.seqinfo is not None

    gr.seqinfo = seq_obj
    assert gr.seqinfo is not None


def test_gr_method_trim():
    gr.seqinfo = seq_obj

    trimmed_gr = gr.trim()

    assert trimmed_gr is not None
    assert (trimmed_gr.start == np.array([101, 102, 103, 104, 105])).all()
    assert (trimmed_gr.width == np.array([11, 21, 16, 30, 5])).all()
