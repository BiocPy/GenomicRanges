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
    "seqnames": [
        "chr1",
        "chr2",
        "chr3",
    ],
    "seqlengths": [110, 112, 118],
    "is_circular": [True, True, False],
    "genome": "hg19",
}


def test_gr_seqInfo():
    gr = genomicranges.from_pandas(df_gr)
    assert gr is not None
    assert gr.seq_info is None

    gr.seq_info = SeqInfo(seq_obj)

    assert gr.seq_info is not None


def test_gr_method_trim():
    gr = genomicranges.from_pandas(df_gr)
    gr.seq_info = SeqInfo(seq_obj)

    trimmed_gr = gr.trim()

    assert trimmed_gr is not None
    assert trimmed_gr.dims == (9, 6)
