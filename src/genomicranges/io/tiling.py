from typing import MutableMapping, Optional, Union
import math
import pandas as pd

from ..utils import split_intervals
from ..GenomicRanges import GenomicRanges
from ..SeqInfo import SeqInfo
from .pdf import fromPandas

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def tileGenome(
    seqlengths: Union[MutableMapping, SeqInfo],
    n: Optional[int] = None,
    width: Optional[int] = None,
) -> "GenomicRanges":
    """Create new genomic regions by partitioning a specified genome.

    Args:
        seqlengths (Union[MutableMapping, SeqInfo]): sequence lengths of each chromosome.
        n (int, optional): number of intervals to split into. Defaults to None.
        width (int, optional): width of each interval. Defaults to None.

    Raises:
        ValueError: either `n` or `width` must be provided but not both.

    Returns:
        GenomicRanges: a new `GenomicRanges` with the tile regions.
    """
    if n is not None and width is not None:
        raise ValueError("either n or width must be provided but not both")

    seqlen_ = seqlengths
    if isinstance(seqlengths, SeqInfo):
        seqlen_ = seqlengths.seqlengths

    all_intervals = []
    for chrm, chrlen in seqlen_.items():
        twidth = None
        if n is not None:
            twidth = math.ceil(chrlen / n)
        elif width is not None:
            twidth = width

        all_intervals.extend(split_intervals(chrm, "*", 1, chrlen, twidth))

    columns = ["seqnames", "strand", "starts", "ends"]
    final_df = pd.DataFrame.from_records(all_intervals, columns=columns)
    final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
    return fromPandas(final_df)
