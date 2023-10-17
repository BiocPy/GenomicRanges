from typing import Dict, List, Optional, Union

from biocframe import BiocFrame
from biocgenerics.combine import combine
from biocgenerics.combine_cols import combine_cols
from biocgenerics.combine_rows import combine_rows
from biocutils import is_list_of_type

from .GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRangesList:
    """Just as it sounds, a `GenomicRangesList` is a named-list like object.

    If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple
    genomic elements, usually where the genes start and end. Genes are themselves made of many sub
    regions, e.g. exons. `GenomicRangesList` allows us to represent this nested structure.

    Currently, this class is limited in functionality, purely a read-only class with basic accessors.

    Typical usage:

    To construct a **GenomicRangesList** object, simply pass in a list of
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` objects and Optionally ``names``.

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

        grl = GenomicRangesList(ranges=[gr1, gr2], names=["gene1", "gene2"])

    Additionally, you may also provide metadata about the genomic elements in the dictionary
    using mcols attribute.
    """

    def __init__(
        self,
        ranges: Union[GenomicRanges, List[GenomicRanges]],
        range_lengths: Optional[List[int]] = None,
        names: Optional[List[str]] = None,
        mcols: Optional[BiocFrame] = None,
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
            mcols (BiocFrame, optional): Metadata about each genomic element. Defaults to None.
            metadata (Optional[Dict], optional): Additional metadata. Defaults to None.
        """
        self._validate(ranges)
        self._data = {"ranges": ranges}

        if range_lengths is None:
            range_lengths = [len(x) for x in ranges]

        self._range_lengths = range_lengths

        if mcols is None:
            mcols = BiocFrame(number_of_rows=len(range_lengths))
        else:
            if not isinstance(mcols, BiocFrame):
                raise TypeError("`mcols` must be a BiocFrame object.")

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

    def __repr__(self):
        from io import StringIO

        from rich.console import Console
        from rich.table import Table

        table = Table(
            title=f"GenomicRangesList with {len(self)} genomic elements",
            show_header=True,
            box=None,
            title_justify="left",
        )

        _rows = []
        rows_to_show = 2
        _top = len(self)
        if _top > rows_to_show:
            _top = rows_to_show

        # top two rows
        for r in range(_top):
            _elem = r
            if self.names is not None:
                _elem = self.names[r]

            _rows.append(f"Genomic element: [bold]{_elem}")
            _rows.append(str(self.ranges[r]))

        if len(self) > rows_to_show:
            if len(self) > 2 * rows_to_show:
                # add ...
                _rows.append("...")

            _last = len(self) - _top
            if _last <= rows_to_show:
                _last = len(self) - _top

            # last set of rows
            for r in range(_last + 1, len(self)):
                _elem = r
                if self.names is not None:
                    _elem = self.names[r]

                _rows.append(f"Genomic element: [bold]{_elem}")
                _rows.append(str(self.ranges[r]))

        for _row in _rows:
            table.add_row(_row)

        console = Console(file=StringIO())
        with console.capture() as capture:
            console.print(table)

        return capture.get()

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
    def mcols(self) -> Optional[BiocFrame]:
        """Get metadata across all genomic elements.

        Returns:
            (BiocFrame, optional): Metadata :py:class:`~biocframe.BiocFrame.Biocframe` or
            None if there is no element level metadata.
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

    def to_pandas(self) -> "DataFrame":
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

    def __getitem__(
        self, args: Union[str, int, tuple, list, slice]
    ) -> Union[GenomicRanges, "GenomicRangesList"]:
        """Access individual genomic elements.

        Args:
            args (Union[str, int, tuple, list, slice]): Name of the genomic element to access.

                Alternatively, if names of genomic elements are not available, you may
                provide an index position of the genomic element to access.

                Alternatively, ``args`` may also specify a list of positions to slice specified either as a
                :py:class:`~list` or :py:class:`~slice`.

                A tuple may also be specified along each dimension. Currently if the tuple contains more than
                one dimension, its ignored.

        Raises:
            TypeError: If ``args`` is not a supported slice argument.

        Returns:
            Union[GenomicRanges, GenomicRangesList]: The genomic element or a new GenomicRangesList of the slice.
        """
        if isinstance(args, int):
            return self.ranges[args]
        elif isinstance(args, str):
            if self.names is not None:
                _idx = self.names.index(args)
                return self.ranges[_idx]
        else:
            new_ranges = None
            new_range_lengths = None
            new_names = None
            new_mcols = None
            new_metadata = self.metadata

            if isinstance(args, tuple):
                # TODO: should figure out what to do with the second dimension later.
                if len(args) >= 1:
                    args = args[0]

            if isinstance(args, list):
                if is_list_of_type(args, bool):
                    if len(args) != len(self):
                        raise ValueError(
                            "`indices` is a boolean vector, length should match the size of the data."
                        )

                    args = [i for i in range(len(args)) if args[i] is True]

                new_ranges = [self.ranges[i] for i in args]
                new_range_lengths = [self._range_lengths[i] for i in args]
                if self.names is not None:
                    new_names = [self.names[i] for i in args]

                if self.mcols is not None:
                    new_mcols = self.mcols[args, :]
            elif isinstance(args, slice):
                new_ranges = self.ranges[args]
                new_range_lengths = self._range_lengths[args]
                if self.names is not None:
                    new_names = self.names[args]

                if self.mcols is not None:
                    new_mcols = self.mcols[args, :]
            else:
                raise TypeError("Arguments to slice is not a list of supported types.")

            return GenomicRangesList(
                new_ranges, new_range_lengths, new_names, new_mcols, new_metadata
            )

        raise TypeError("Arguments to slice is not supported.")

    def __len__(self) -> int:
        """Number of genomic elements.

        Returns:
            int: number of genomic elements.
        """
        return len(self._range_lengths)

    @classmethod
    def empty(cls, n: int):
        """Create an empty ``n``-length `GenomicRangesList` object.

        Returns:
            same type as caller, in this case a `GenomicRangesList`.
        """
        _range_lengths = [0] * n

        return cls(ranges=GenomicRanges.empty(), range_lengths=_range_lengths)

    def combine(self, *other: "GenomicRangesList") -> "GenomicRangesList":
        """Combine multiple `GenomicRangesList` objects by row.

        Note: Fills missing columns with an array of `None`s.

        Args:
            *other (GenomicRangesList): Objects to combine.

        Raises:
            TypeError: If all objects are not of type GenomicRangesList.

        Returns:
            GenomicRangesList: A combined `GenomicRangesList` object.
        """
        if not is_list_of_type(other, GenomicRangesList):
            raise TypeError("All objects to combine must be GenomicRangesList objects.")

        # ranges: Union[GenomicRanges, List[GenomicRanges]],
        # range_lengths: Optional[List[int]] = None,
        # names: Optional[List[str]] = None,
        # mcols: Optional[BiocFrame] = None,
        # metadata: Optional[dict] = None,

        all_objects = [self] + list(other)

        new_ranges = combine(*[obj.ranges for obj in all_objects])
        new_names = combine(*[obj.names for obj in all_objects])
        new_mcols = combine(*[obj.mcols for obj in all_objects])
        new_metadata = {}
        for i, obj in enumerate(all_objects):
            if obj.metadata is not None:
                new_metadata[i] = obj.metadata

        current_class_const = type(self)
        return current_class_const(
            new_ranges, names=new_names, mcols=new_mcols, metadata=new_metadata
        )


@combine.register(GenomicRangesList)
def _combine_grl(*x: GenomicRangesList):
    if not is_list_of_type(x, GenomicRangesList):
        raise ValueError(
            "All elements to `combine` must be `GenomicRangesList` objects."
        )

    return x[0].combine(*x[1:])


@combine_rows.register(GenomicRangesList)
def _combine_rows_grl(*x: GenomicRangesList):
    if not is_list_of_type(x, GenomicRangesList):
        raise ValueError(
            "All elements to `combine_rows` must be `GenomicRangesList` objects."
        )

    return x[0].combine(*x[1:])


@combine_cols.register(GenomicRangesList)
def _combine_cols_grl(*x: GenomicRangesList):
    if not is_list_of_type(x, GenomicRangesList):
        raise ValueError(
            "All elements to `combine_cols` must be `GenomicRangesList` objects."
        )

    raise NotImplementedError(
        "`combine_cols` is not implemented for `GenomicRangesList` objects."
    )
