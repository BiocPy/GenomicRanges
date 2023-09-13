from collections import UserDict
from typing import Dict, List, Mapping, MutableMapping, Optional, Union

from biocframe import BiocFrame
from pandas import DataFrame, concat

from .GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

BiocOrPandasFrame = Union[DataFrame, BiocFrame]


class GenomicRangesList(UserDict):
    """Just as it sounds, a ``GenomicRangesList`` is a dictionary, where the keys represent features and the value a
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` object.

    If you are wondering why you need this class, a :py:class:`genomicranges.GenomicRanges.GenomicRanges`
    object lets us specify mutiple genomic elements, usually where the genes start and end. Genes are themselves made
    of many sub regions, e.g. exons. ``GenomicRangesList`` allows us to represent this nested structure.

    Currently this class is limited in functionality. Purely a read-only class with basic
    accessors.

    Typical usage example:

    To construct a **GenomicRangesList** object, simply pass in a dictionary

    .. code-block:: python

        gr1 = GenomicRanges(
            {
                "seqnames": ["chr1", "chr2", "chr1", "chr3"],
                "starts": [1, 3, 2, 4],
                "ends": [10, 30, 50, 60],
                "strand": ["-", "+", "*", "+"],
                "score": [1, 2, 3, 4],
            }
        )

        gr2 = GenomicRanges(
            {
                "seqnames": ["chr2", "chr4", "chr5"],
                "starts": [3, 6, 4],
                "ends": [30, 50, 60],
                "strand": ["-", "+", "*"],
                "score": [2, 3, 4],
            }
        )

        grl = GenomicRangesList(gene1=gr1, gene2=gr2)

    Additionally you may also provide metadata about the genomic elements in the dictionary
    using :py:attr:`~genomicranges.GenomicRangesList.GenomicRangesList.mcols`.
    """

    def __init__(
        self,
        mcols: BiocOrPandasFrame = None,
        metadata: Optional[Mapping] = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)

        if mcols is None:
            mcols = BiocFrame(number_of_rows=len(self), row_names=list(self.keys()))

        self._validate_mcols(mcols)
        self._mcols = mcols

        self._metadata = metadata

    def _validate(self):
        """Internal wrapper method to validate the object."""
        # validate assays to make sure they are have same dimensions
        self._validate_mcols(self._mcols)

    def _validate_mcols(self, mcols):
        """Internal method to validate genomic elements (mcols)."""
        if len(self) != len(mcols):
            raise ValueError("Number of elements and mcols are not the same length!")

    @property
    def metadata(self) -> Optional[MutableMapping]:
        """Get metadata.

        Returns:
            Optional[MutableMapping]: Metadata object, usually a dictionary.
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional[MutableMapping]):
        """Set metadata.

        Args:
            metadata (Optional[MutableMapping]): New metadata object.
        """
        self._metadata = metadata

    @property
    def mcols(self) -> BiocOrPandasFrame:
        """Get metadata across all genomic elements.

        Returns:
            BiocOrPandasFrame: Metadata frame.
        """

        return self._mcols

    def __setitem__(self, key, value):
        if not isinstance(value, GenomicRanges):
            raise TypeError("Value must be of type `GenomicRanges`.")

        super().__setitem__(key, value)

    def element_nrows(self) -> List[int]:
        """Get a vector of the length of each element.

        Returns:
            List[int]: An integer vector where each value corresponds to the length of
            the contained GenomicRanges object.
        """
        return self._generic_accessor("__len__", func=True)

    def is_empty(self) -> bool:
        """Whether ``GRangesList`` has no elements or if all its elements are empty.

        Returns:
            bool: True if the object has no elements.
        """
        if len(self) == 0:
            return True

        element_lengths = self.element_nrows()
        if all([True if x == 0 else False for x in element_lengths]):
            return True

        return False

    def _generic_accessor(self, prop: str, func: bool = False) -> Dict[str, List]:
        _all_prop = {}
        for k, v in self.items():
            _val = getattr(v, prop)

            if func is True:
                _val = _val()

            _all_prop[k] = _val

        return _all_prop

    # TODO: convert some of these properties to a factorized array

    @property
    def seqnames(self) -> Dict[str, List[str]]:
        """Get all sequence names.

        Returns:
            Dict[str, List[str]]: A list with the same length as keys in the object,
            each element in the list contains another list of sequence names.
        """
        return self._generic_accessor("seqnames")

    def ranges(self) -> Dict[str, List[str]]:
        """Get all ranges.

        Returns:
            Dict[str, List[str]]: A list with the same length as keys in the object,
            each element in the list contains another list of ranges names.
        """
        return self._generic_accessor("ranges", func=True)

    @property
    def start(self) -> Dict[str, List[int]]:
        """Get all start positions.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("start")

    @property
    def end(self) -> Dict[str, List[int]]:
        """Get all end positions.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("end")

    @property
    def width(self) -> Dict[str, List[int]]:
        """Get width of all regions across all elements.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("width")

    @property
    def strand(self) -> Dict[str, List[int]]:
        """Get strand of all regions across all elements.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("strand")

    @property
    def seq_info(self) -> Dict[str, List[int]]:
        """Get information about the underlying sequences.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("seq_info")

    @property
    def is_circular(self) -> Dict[str, List[int]]:
        """Get the circularity flag.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("is_circular")

    @property
    def genome(self) -> Dict[str, List[int]]:
        """Get genome of the underlying sequences.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("genome")

    @property
    def score(self) -> Dict[str, List[int]]:
        """Get score about the underlying sequences.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("score")

    def to_pandas(self) -> DataFrame:
        """Coerce object to a :py:class:`pandas.DataFrame`.

        Returns:
            DataFrame: A :py:class:`~pandas.DataFrame` object.
        """
        all_index = []
        all_dfs = []

        for k, v in self.items():
            _idx = [k] * len(v)
            all_index.extend(_idx)
            all_dfs.append(v.to_pandas())

        all_concat = concat(all_dfs)
        all_concat.index = all_index

        return all_concat

    def add_element(self, key, value, element_metadata):
        raise NotImplementedError("Adding new elements is not yet implemented!")
