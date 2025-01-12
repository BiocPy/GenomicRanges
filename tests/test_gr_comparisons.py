from random import random

import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges.GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

gr = GenomicRanges(
    seqnames=["chr1", "chr2", "chr2", "chr2", "chr1", "chr1", "chr3", "chr3", "chr3", "chr3"],
    ranges=IRanges([x for x in range(1, 11)], [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(1, 11),
            "GC": [random() for _ in range(10)],
        }
    ),
)


def test_matches():
    assert gr is not None

    query_hits = gr.match(gr[2:5])

    assert query_hits is not None
    assert np.all(query_hits == [2, 3, 4])


def test_order():
    assert gr is not None

    order = gr.order()
    assert (order == np.array([5, 0, 4, 1, 2, 3, 6, 7, 8, 9])).all()


def test_sort():
    assert gr is not None
    result = gr.sort()

    assert result is not None

    result = gr.sort(decreasing=True)
    assert result is not None


def test_rank():
    assert gr is not None
    result = gr.rank()

    assert result is not None
    assert (result == np.array([1, 3, 4, 5, 2, 0, 6, 7, 8, 9])).all()
