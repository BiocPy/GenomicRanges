from random import random

import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from genomicranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

subject = GenomicRanges(
    seqnames=[
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
    ranges=IRanges(range(1, 11), [10] * 10),
    strand=["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame(
        {
            "score": range(0, 10),
            "GC": [random() for _ in range(10)],
        }
    ),
)

query = GenomicRanges(
    seqnames=["chr2", "chr2"],
    ranges=IRanges([4, 3], [3, 4]),
    strand=["+", "+"],
)


def test_find_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query)

    assert res is not None
    assert isinstance(res, BiocFrame)
    assert np.all(res.get_column("self_hits") == [1, 2, 3, 1, 2, 3])
    assert np.all(res.get_column("query_hits") == [0, 0, 0, 1, 1, 1])

    res = query.find_overlaps(subject, ignore_strand=True)

    assert res is not None
    assert isinstance(res, BiocFrame)
    assert np.all(res.get_column("query_hits") == [1, 1, 2, 2, 3, 3])
    assert np.all(res.get_column("self_hits") == [0, 1, 0, 1, 0, 1])


def test_find_overlaps_query_type():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query, query_type="within")

    assert res is not None
    assert np.all(res.get_column("self_hits") == [1, 2, 3, 1, 2])
    assert np.all(res.get_column("query_hits") == [0, 0, 0, 1, 1])


def test_find_overlaps_threads():
    assert subject is not None
    assert query is not None

    res = subject.find_overlaps(query, query_type="within", num_threads=3)

    assert res is not None
    assert np.all(res.get_column("self_hits") == [1, 2, 3, 1, 2])
    assert np.all(res.get_column("query_hits") == [0, 0, 0, 1, 1])


def test_find_overlaps_rtrip():
    x = GenomicRanges(["chr1", "chr1"], IRanges([2, 9], [7, 19]), strand=["+", "-"])
    y = GenomicRanges(["chr1"], IRanges([5], [10]), strand=["*"])

    resxy = x.find_overlaps(y)
    assert resxy is not None
    assert np.all(resxy.get_column("self_hits") == [0, 1])
    assert np.all(resxy.get_column("query_hits") == [0, 0])

    resyx = y.find_overlaps(x)
    assert resyx is not None
    assert np.all(resyx.get_column("self_hits") == [0, 0])
    assert np.all(resyx.get_column("query_hits") == [0, 1])


def test_count_overlaps():
    assert subject is not None
    assert query is not None

    res = subject.count_overlaps(query)

    assert res is not None
    assert isinstance(res, np.ndarray)
    assert np.all(res == [3, 3])


def test_subset_by_overlaps():
    query2 = GenomicRanges(
        seqnames=["chr2", "chr4"],
        ranges=IRanges(start=[4, 4], width=[2, 2]),
        strand=["+", "-"],
    )

    assert subject is not None
    assert query2 is not None

    res = subject.subset_by_overlaps(query2)

    assert res is not None
    assert isinstance(res, GenomicRanges)
    assert len(res) == 3
