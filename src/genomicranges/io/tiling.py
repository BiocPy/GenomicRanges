import math
from typing import Dict, Optional, Union

import pandas as pd

from ..interval import split_intervals

# from ..GenomicRanges import GenomicRanges
from ..Seqinfo import Seqinfo
from .pdf import from_pandas

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def tile_genome(
    seqlengths: Union[Dict, Seqinfo],
    n: Optional[int] = None,
    width: Optional[int] = None,
) -> "GenomicRanges":
    """Create new genomic regions by partitioning a specified genome.

    If ``n`` is provided, the region is split into ``n`` intervals. The last interval may
    not contain the same 'width' as the other regions.

    Alternatively, ``width`` may be provided for each interval. Similarly, the last region
    may be less than ``width``.

    Either ``n`` or ``width`` must be provided but not both.

    Args:
        seqlengths (Union[Dict, Seqinfo]): Sequence lengths of each chromosome.

            ``seqlengths`` may be a dictionary, where keys specify the chromosome and the value is
            thelength of each chromosome in the genome.

            Alternatively, ``seqlengths`` may be an instance of
            :py:class:`~genomicranges.Seqinfo.Seqinfo`.

        n (int, optional): Number of intervals to split into.
            Defaults to None, then 'width' of each interval is computed from ``seqlengths``.

        width (int, optional): Width of each interval. Defaults to None.

    Raises:
        ValueError: Either ``n`` or ``width`` must be provided but not both.

    Returns:
        GenomicRanges: The genome with the tiled regions.
    """

    if n is not None and width is not None:
        raise ValueError("Both `n` or `width` are provided!")

    seqlen_ = seqlengths
    if isinstance(seqlengths, Seqinfo):
        seqlen_ = seqlengths.seqlengths()

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
    return from_pandas(final_df)
