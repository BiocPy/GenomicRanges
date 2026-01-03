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

q_gr = GenomicRanges(
    seqnames=["chr2", "chr2"],
    ranges=IRanges([4, 3], [3, 4]),
    strand=["+", "+"],
)


def test_nearest():
    assert gr is not None

    query_hits = gr.nearest(q_gr)

    assert query_hits is not None
    assert np.all(query_hits == [1, 1])  # R returns [3,3], select is arbitrary so its ok

    query_hits = q_gr.nearest(gr)
    assert np.all(
        query_hits == [None, 1, 1, 1, None, None, None, None, None, None]
    )  # R has 1,1,1 or 0,0,0 but both are fine since they are near

    query_hits = gr.nearest(q_gr, select="all")
    assert np.all(query_hits.get_column("query_hits") == [0, 0, 0, 1, 1, 1])
    assert np.all(query_hits.get_column("self_hits") == [1, 2, 3, 1, 2, 3])

    query_hits = gr.nearest(q_gr, ignore_strand=True)
    assert np.all(query_hits == [1, 1])  # R returns [3,3], select is arbitrary so its ok

    query_hits = gr.nearest(q_gr, select="all", ignore_strand=True)
    assert np.all(query_hits.get_column("query_hits") == [0, 0, 0, 1, 1, 1])
    assert np.all(query_hits.get_column("self_hits") == [1, 2, 3, 1, 2, 3])

    query_hits = gr.nearest(q_gr, select="all", ignore_strand=True, num_threads=3)
    assert np.all(query_hits.get_column("query_hits") == [0, 0, 0, 1, 1, 1])
    assert np.all(query_hits.get_column("self_hits") == [1, 2, 3, 1, 2, 3])


def test_precede():
    assert gr is not None

    query_hits = gr.precede(q_gr)

    assert query_hits is not None
    assert np.all(query_hits == [None, None])

    query_hits = q_gr.precede(gr)
    assert np.all(query_hits == [None] * len(gr))

    query_hits = gr.precede(q_gr, select="all")
    assert np.all(query_hits.get_column("query_hits") == [])
    assert np.all(query_hits.get_column("self_hits") == [])

    query_hits = gr.precede(q_gr, select="all", num_threads=3)
    assert np.all(query_hits.get_column("query_hits") == [])
    assert np.all(query_hits.get_column("self_hits") == [])


def test_follow():
    assert gr is not None

    query_hits = gr.follow(q_gr)

    assert query_hits is not None
    assert np.all(query_hits == [None, None])

    query_hits = q_gr.follow(gr)
    assert np.all(query_hits == [None] * len(gr))

    query_hits = gr.follow(q_gr, select="all")
    assert np.all(query_hits.get_column("query_hits") == [])
    assert np.all(query_hits.get_column("self_hits") == [])

    query_hits = gr.follow(q_gr, select="all", num_threads=3)
    assert np.all(query_hits.get_column("query_hits") == [])
    assert np.all(query_hits.get_column("self_hits") == [])


def test_distance():
    query = GenomicRanges(["A", "A", "A"], ranges=IRanges([1, 2, 10], [5, 7, 2]))
    subject = GenomicRanges(["A", "A", "A"], ranges=IRanges([6, 5, 13], [5, 6, 3]))

    distance = subject.distance(query)

    assert distance is not None
    assert isinstance(distance, np.ndarray)
    assert distance.tolist() == [0, 0, 1]
