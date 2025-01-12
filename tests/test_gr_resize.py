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


def test_resize_default():
    assert gr is not None

    resized_gr = gr.resize(width=10)

    assert resized_gr is not None
    assert isinstance(resized_gr, GenomicRanges)
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1, 2, 3, 4, 5, 6, 7, 8, 1, 1]))
    assert np.all(resized_gr.end == np.array([10, 11, 12, 13, 14, 15, 16, 17, 10, 10]))


def test_resize_fix_end():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end")

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1, 1, 1, 1, 1, 1, 1, 1, 9, 10]))
    assert np.all(resized_gr.end == np.array([10, 10, 10, 10, 10, 10, 10, 10, 18, 19]))


def test_resize_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, ignore_strand=True)

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
    assert np.all(resized_gr.end == np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19]))


def test_resize_fix_end_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end", ignore_strand=True)

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1] * 10))
    assert np.all(resized_gr.end == np.array([10] * 10))


def test_resize_fix_center_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center", ignore_strand=True)

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5]))
    assert np.all(resized_gr.end == np.array([10, 10, 11, 11, 12, 12, 13, 13, 14, 14]))


def test_resize_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center")

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5]))
    assert np.all(resized_gr.end == np.array([10, 10, 11, 11, 12, 12, 13, 13, 14, 14]))


def test_resize_default_width_odd():
    assert gr is not None

    resized_gr = gr.resize(width=11)

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([0, 2, 3, 4, 5, 6, 7, 8, 0, 0]))
    assert np.all(resized_gr.end == np.array([10, 12, 13, 14, 15, 16, 17, 18, 10, 10]))


def test_resize_default_width_odd_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=11, fix="center")

    assert resized_gr is not None
    assert resized_gr.get_strand(as_type="list") == ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"]
    assert np.all(resized_gr.start == np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5]))
    assert np.all(resized_gr.end == np.array([10, 11, 11, 12, 12, 13, 13, 14, 14, 15]))
