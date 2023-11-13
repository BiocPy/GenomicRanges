from typing import (
    List,
    Optional,
    Union,
)

import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from .SeqInfo import SeqInfo

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges:
    """`GenomicRanges` provides a container class to represent and operate over genomic regions and annotations.

    Additionally, `GenomicRanges` may also contain `Sequence Information` (checkout
    :py:class:`~genomicranges.SeqInfo.SeqInfo`) as part of its metadata. It contains for each
    sequence name (or chromosome) in the gene model, its length. Additionally, (checkout
    :py:class:`~genomicranges.SeqInfo.SeqInfo`) might also contain metadata about the
    genome, e.g. if it's circular (`is_circular`) or not.

    Note: The documentation for some of the methods are derived from the
    `GenomicRanges R/Bioconductor package <https://github.com/Bioconductor/GenomicRanges>`_.
    """

    required_columns = ["seqnames", "starts", "ends", "strand"]

    def __init__(
        self,
        seqnames: Union[List[str], np.ndarray],
        ranges: IRanges,
        strand: Optional[List[str, np.ndarray]] = None,
        names: Optional[List] = None,
        mcols: Optional[BiocFrame] = None,
        seq_info: Optional[SeqInfo] = None,
        metadata: Optional[dict] = None,
        validate: bool = True,
    ):
        pass
