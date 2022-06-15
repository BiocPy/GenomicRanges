import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_ucsc():
    gr = GenomicRanges.fromUCSC("hg19")
    assert gr is not None
    assert gr.len() > 0
    assert len(gr) == gr.len()
    assert gr.mcols() is not None
    assert gr.mcols().shape[0] > 0
    assert gr.granges() is not None
