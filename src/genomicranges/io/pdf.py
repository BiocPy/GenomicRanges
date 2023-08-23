from collections import OrderedDict

from pandas import DataFrame

# from ..GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def from_pandas(data: DataFrame) -> "GenomicRanges":
    """Read a :py:class:`~pandas.DataFrame` into :py:class:`~genomicranges.GenomicRanges.GenomicRanges`.

    Args:
        data (DataFrame): `DataFrame` object with genomic positions.
            Must contain 'seqnames', 'starts' & 'ends' columns.

    Returns:
        GenomicRanges: object.
    """

    from ..GenomicRanges import GenomicRanges

    obj = OrderedDict()

    for col in data.columns:
        obj[col] = data[col].to_list()

    rindex = None
    if data.index is not None:
        rindex = data.index.to_list()

    return GenomicRanges(obj, row_names=rindex, column_names=data.columns.to_list())
