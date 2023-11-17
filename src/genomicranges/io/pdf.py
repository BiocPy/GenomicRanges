from warnings import warn

# from ..GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def from_pandas(data: "pandas.DataFrame") -> "GenomicRanges":
    """Alias to :py:meth:`genomicranges.GenomicRanges.GenomicRanges.from_pandas`.

    Returns:
        A ``GenomicRanges`` object representing intervals.
    """

    warn(
        "This method is deprecated. Use 'GenomicRanges.from_pandas' instead.",
        UserWarning,
    )

    from ..GenomicRanges import GenomicRanges

    return GenomicRanges.from_pandas(data)
