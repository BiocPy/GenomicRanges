import pytest
import pandas as pd
from genomicranges.GenomicRanges import GenomicRanges
from random import random

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

@pytest.mark.skip(reason="takes too long")
def test_ucsc():
    gr = GenomicRanges.fromUCSC("hg19")
    assert gr is not None
    assert len(gr) > 0
    assert gr.granges() is not None
