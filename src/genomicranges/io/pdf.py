import pandas as pd
from collections import OrderedDict

from ..GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def fromPandas(data: pd.DataFrame) -> GenomicRanges:
    """Convert a pandas `Dataframe` to `GenomicRanges`.

    Args:
        data (pd.DataFrame): a Pandas `DataFrame` object containing genomic positions.
            Must contain `seqnames`, `starts` & `ends` columns.

    Returns:
        GenomicRanges: A `GenomicRanges` object representing genomic positions.
    """

    obj = OrderedDict()

    for col in data.columns:
        obj[col] = data[col].to_list()

    rindex = None
    if data.index is not None:
        rindex = data.index.to_list()

    return GenomicRanges(obj, rowNames=rindex, columnNames=data.columns.to_list())
