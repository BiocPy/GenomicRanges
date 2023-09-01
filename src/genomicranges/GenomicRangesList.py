from collections import UserDict
from typing import Dict, List

from .GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRangesList(UserDict):
    """Just as it sounds, a ``GenomicRangesList`` is a dictionary, where the keys represent features and the value a
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` object.

    If you are wondering why I need this class, a :py:class:`genomicranges.GenomicRanges.GenomicRanges`
    object lets us specify mutiple regions, usually where the genes start and end. Genes are themselves made
    of many sub regions, e.g. exons. ``GenomicRangesList`` allows us to represent this nested structure.

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
    """

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
            print(k, v)
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

    def mcols(self) -> Dict[str, List[int]]:
        """Get all metadata columns.

        Returns:
            Dict[str, List[int]]: A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("mcols", func=True)

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

    def as_data_frame(self):
        """Get the object as data frame.

        Raises:
            NotImplementedError: Not yet implemented.
        """
        raise NotImplementedError("coercion to dataframe is not available yet.")
