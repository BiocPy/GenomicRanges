import pandas as pd
from genomicranges.SeqInfo import SeqInfo
from random import random
import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


df_gr = pd.DataFrame(
    {
        "seqnames": [
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
        "starts": range(100, 110),
        "ends": range(110, 120),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

seq_obj = {
    "seqnames": ["chr1", "chr2", "chr3",],
    "seqlengths": [110, 112, 118],
    "isCircular": [True, True, False],
    "genome": "hg19",
}


def test_gr_seqInfo():
    gr = genomicranges.fromPandas(df_gr)
    assert gr is not None
    assert gr.seqInfo is None

    gr.seqInfo = SeqInfo(seq_obj)

    assert gr.seqInfo is not None


def test_gr_method_trim():
    gr = genomicranges.fromPandas(df_gr)
    gr.seqInfo = SeqInfo(seq_obj)

    trimmed_gr = gr.trim()

    assert trimmed_gr is not None
    assert trimmed_gr.dims == (9, 6)
