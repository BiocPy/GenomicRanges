import pytest
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random
import pandas as pd
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_create_gr():
    gr = GenomicRanges(
        seqnames=["chr1"] * 10,
        ranges=IRanges(start=range(100, 110), width=range(110, 120)),
    )

    assert gr is not None
    assert gr._seqnames.dtype == np.uint8

    gr16 = GenomicRanges(
        seqnames=[f"chr{i}" for i in range(500)],
        ranges=IRanges(start=range(0, 500), width=range(10, 510)),
    )

    assert gr16 is not None
    assert gr16._seqnames.dtype == np.uint16

    gr32 = GenomicRanges(
        seqnames=[f"chr{i}" for i in range(2**16 + 1)],
        ranges=IRanges(start=range(0, 2**16 + 1), width=range(10, 2**16 + 11)),
    )

    assert gr32 is not None
    assert gr32._seqnames.dtype == np.uint32
