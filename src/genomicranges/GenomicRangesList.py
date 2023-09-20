from typing import Dict, List, Optional, Union

from biocframe import BiocFrame
from pandas import DataFrame

from .GenomicRanges import GenomicRanges
from .utils import is_list_of_type

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

BiocOrPandasFrame = Union[DataFrame, BiocFrame]


class GenomicRangesList:
    """Just as it sounds, a `GenomicRangesList` is a named-list like object.

    If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple
    genomic elements, usually where the genes start and end. Genes are themselves made of many sub
    regions, e.g. exons. `GenomicRangesList` allows us to represent this nested structure.

    Currently, this class is limited in functionality, purely a read-only class with basic accessors.

    Typical usage example:

    To construct a `GenomicRangesList` object, simply pass in a list of
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` objects and Optionally ``names``.

    ```python
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

        grl = GenomicRangesList(ranges=[gr1, gr2], names=["gene1", "gene2"])
    ```

    Additionally, you may also provide metadata about the genomic elements in the dictionary
    using mcols attribute.
    """

    def __init__(
        self,
        ranges: Union[GenomicRanges, List[GenomicRanges]],
        range_lengths: Optional[List[int]] = None,
        names: Optional[List[str]] = None,
        mcols: BiocOrPandasFrame = None,
        metadata: Optional[dict] = None,
    ):
        """Initialize a `GenomicRangesList` object.

        Args:
            ranges (GenomicRangesList, optional): List of genomic elements.
                All elements in this list must be :py:class:`genomicranges.GenomicRanges.GenomicRanges`
                class.
            range_lengths (Optional[List[int]], optional): Number of ranges within each genomic element.
                Defaults to None.
            names (Optional[List[str]], optional): Names of the genomic elements.
                The length of this must match the number of genomic elements in ``ranges``.
                Defaults to None.
            mcols (BiocOrPandasFrame, optional): Metadata about each genomic element. Defaults to None.
            metadata (Optional[Dict], optional): Additional metadata. Defaults to None.
        """
        self._validate(ranges)
        self._data = {"ranges": ranges}

        if range_lengths is None:
            range_lengths = [len(x) for x in ranges]

        self._range_lengths = range_lengths

        if mcols is None:
            mcols = BiocFrame(number_of_rows=len(range_lengths))

        self._data["mcols"] = mcols

        self._names = names
        self._metadata = {} if metadata is None else metadata

    def _validate(self, ranges: Union[GenomicRanges, List[GenomicRanges]]):
        if not (
            isinstance(ranges, GenomicRanges) or is_list_of_type(ranges, GenomicRanges)
        ):
            raise TypeError(
                "`ranges` must be either a `GenomicRanges` or a list of `GenomicRanges`."
            )

    @property
    def metadata(self) -> dict:
        """Get metadata.

        Returns:
            dict: Metadata object, usually a dictionary.
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional[dict]):
        """Set metadata.

        Args:
            metadata (dict, optional): New metadata object.
        """
        self._metadata = {} if metadata is None else metadata

    @property
    def ranges(self) -> Dict[str, List[str]]:
        """Get all ranges.

        Returns:
            Dict[str, List[str]]: A list with the same length as keys in the object,
            each element in the list contains another list of ranges names.
        """
        return self._data["ranges"]

    @property
    def groups(self) -> Optional[list]:
        """Get names of all genomic elements.

        Returns:
            (list, optional): A list with the same length as
            :py:attr:`genomicranges.GenomicRanges.GenomicRangesList.ranges`,
            with each element specifying the name of the element. None if names are not provided.
        """
        return self._names

    @property
    def names(self) -> Optional[list]:
        """Alias to :py:attr:`genomicranges.GenomicRanges.GenomicRangesList.groups`.

        Returns:
            (list, optional): A list with the same length as
            :py:attr:`genomicranges.GenomicRanges.GenomicRangesList.ranges`,
            with each element specifying the name of the element. None if names are not provided.
        """
        return self.groups

    @property
    def mcols(self) -> Optional[BiocOrPandasFrame]:
        """Get metadata across all genomic elements.

        Returns:
            (BiocOrPandasFrame, optional): Metadata frame or None if there is no element level metadata.
        """
        if "mcols" in self._data:
            return self._data["mcols"]

        return None

    def _generic_accessor(self, prop: str, func: bool = False) -> Dict[str, list]:
        _all_prop = {}
        _ranges = self.ranges
        _groups = self.names

        for i in range(len(_ranges)):
            _val = getattr(_ranges[i], prop)

            if func is True:
                _val = _val()

            _key = i
            if _groups is not None:
                _key = _groups[i]

            _all_prop[_key] = _val

        return _all_prop

    def element_nrows(self) -> Dict[str, List[str]]:
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

    # TODO: convert some of these properties to a factorized array

    @property
    def seqnames(self) -> Dict[str, List[str]]:
        """Get all sequence names.

        Returns:
            Dict[str, List[str]]: A list with the same length as keys in the object,
            each element in the list contains another list of sequence names.
        """
        return self._generic_accessor("seqnames")

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
        from pandas import concat

        all_index = []
        all_dfs = []

        _ranges = self.ranges
        _groups = self.names

        for i in range(len(_ranges)):
            v = _ranges[i]

            _key = i
            if _groups is not None:
                _key = _groups[i]

            _idx = [_key] * len(v)
            all_index.extend(_idx)
            all_dfs.append(v.to_pandas())

        all_concat = concat(all_dfs)
        all_concat.index = all_index

        return all_concat

    def add_element(self, key, value, element_metadata):
        raise NotImplementedError("Adding new elements is not yet implemented!")

    def __getitem__(self, args: Union[str, int]) -> "GenomicRanges":
        """Access individual genomic elements.

        Args:
            args (Union[str, int]): Name of the genomic element to access.

                Alternatively, if names of genomic elements are not available, you may
                provide an index position of the genomic element to access.

        Raises:
            TypeError: If ``args`` is not a string nor integer.

        Returns:
            GenomicRanges: The genomic element.
        """
        if isinstance(args, int):
            return self.ranges[args]
        elif isinstance(args, str):
            if self.names is not None:
                _idx = self.names.index(args)
                return self.ranges[_idx]

        raise TypeError("args must be either a string or an integer.")

    def __len__(self) -> int:
        """Number of genomic elements.

        Returns:
            int: number of genomic elements.
        """
        return len(self._range_lengths)

    @classmethod
    def empty(cls, n: int):
        """Create an ``n``-length `GenomicRangesList` object.

        Returns:
            same type as caller, in this case a `GenomicRangesList`.
        """
        _range_lengths = [0] * n

        return cls(ranges=GenomicRanges.empty(), range_lengths=_range_lengths)
