import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random
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


def test_resize_default():
    assert gr is not None

    resized_gr = gr.resize(width=10)

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 113, 103, 104, 100])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_fix_end():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end")

    assert resized_gr is not None
    assert (resized_gr.start == np.array([102, 102, 118, 124, 105])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, ignore_strand=True)

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 102, 103, 104, 105])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_fix_end_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="end", ignore_strand=True)

    assert resized_gr is not None
    assert (resized_gr.start == np.array([102, 113, 118, 124, 100])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_fix_center_ignore_strand():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center", ignore_strand=True)

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 107, 110, 114, 102])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=10, fix="center")

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 107, 110, 114, 102])).all()
    assert (resized_gr.width == np.array([10] * 5)).all()


def test_resize_default_width_odd():
    assert gr is not None

    resized_gr = gr.resize(width=11)

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 112, 103, 104, 99])).all()
    assert (resized_gr.width == np.array([11] * 5)).all()


def test_resize_default_width_odd_fix_center():
    assert gr is not None

    resized_gr = gr.resize(width=11, fix="center")

    assert resized_gr is not None
    assert (resized_gr.start == np.array([101, 107, 110, 113, 102])).all()
    assert (resized_gr.width == np.array([11] * 5)).all()
