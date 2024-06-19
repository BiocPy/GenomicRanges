from typing import Dict, List, Optional, Sequence, Union
from warnings import warn

import biocutils as ut
from biocframe import BiocFrame
import numpy as np

from .GenomicRanges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _validate_ranges(ranges, num_ranges):
    if ranges is None:
        raise ValueError("'ranges' cannot be None.")

    if not (
        isinstance(ranges, GenomicRanges) or ut.is_list_of_type(ranges, GenomicRanges)
    ):
        raise TypeError(
            "`ranges` must be either a `GenomicRanges` or a list of `GenomicRanges`."
        )

    if isinstance(ranges, list) and sum([len(x) for x in ranges]) != num_ranges:
        raise ValueError(
            "Length of 'ranges' does not match the number of genomic elements.",
            f"Need to be {num_ranges}, provided {len(ranges)}.",
        )
    elif isinstance(ranges, GenomicRanges) and len(ranges) != num_ranges:
        raise ValueError(
            "Length of 'ranges' does not match the number of genomic elements.",
            f"Need to be {num_ranges}, provided {len(ranges)}.",
        )


def _validate_optional_attrs(mcols, names, num_ranges):
    if not isinstance(mcols, BiocFrame):
        raise TypeError("'mcols' is not a `BiocFrame` object.")

    if mcols.shape[0] != num_ranges:
        raise ValueError(
            "Length of 'mcols' does not match the number of genomic elements."
        )

    if names is not None:
        if len(names) != num_ranges:
            raise ValueError(
                "Length of 'names' does not match the number of genomic elements."
            )

        if any(x is None for x in names):
            raise ValueError("'names' cannot contain None values.")


def _sanitize_range_lengths(x):
    if not isinstance(x, np.ndarray):
        return np.array(x, dtype=np.int32)

    return x


class GenomicRangesListIter:
    """An iterator to a :py:class:`~GenomicRangesList` object."""

    def __init__(self, obj: "GenomicRangesList") -> None:
        """Initialize the iterator.

        Args:
            obj:
                source object to iterate.
        """
        self._grl = obj
        self._current_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._current_index < len(self._grl):
            iter_row_index = (
                self._grl.names[self._current_index]
                if self._grl.names is not None
                else None
            )

            iter_slice = self._grl[self._current_index]
            self._current_index += 1
            return (iter_row_index, iter_slice)

        raise StopIteration


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

        a = GenomicRanges(
            seqnames=["chr1", "chr2", "chr1", "chr3"],
            ranges=IRanges([1, 3, 2, 4], [10, 30, 50, 60]),
            strand=["-", "+", "*", "+"],
            mcols=BiocFrame({"score": [1, 2, 3, 4]}),
        )

        b = GenomicRanges(
            seqnames=["chr2", "chr4", "chr5"],
            ranges=IRanges([3, 6, 4], [30, 50, 60]),
            strand=["-", "+", "*"],
            mcols=BiocFrame({"score": [2, 3, 4]}),
        )

        grl = GenomicRangesList(ranges=[gr1, gr2], names=["gene1", "gene2"])

    Additionally, you may also provide metadata about the genomic elements in the dictionary
    using mcols attribute.
    """

    def __init__(
        self,
        ranges: Union[GenomicRanges, List[GenomicRanges]],
        range_lengths: Optional[Sequence[int]] = None,
        names: Optional[List[str]] = None,
        mcols: Optional[BiocFrame] = None,
        metadata: Optional[dict] = None,
        validate: bool = True,
    ):
        """Initialize a `GenomicRangesList` object.

        Args:
            ranges:
                List of genomic elements.
                All elements in this list must be
                :py:class:`genomicranges.GenomicRanges.GenomicRanges` objects.

            range_lengths:
                Number of ranges within each genomic element.
                Defaults to None, and is inferred from ``ranges``.

            names:
                Names of the genomic elements.
                The length of this must match the number of
                genomic elements in ``ranges``. Defaults to None.

            mcols:
                Metadata about each genomic element. Defaults to None.

            metadata:
                Additional metadata. Defaults to None.

            validate:
                Internal use only.
        """
        self._ranges = ranges

        if range_lengths is None:
            if isinstance(ranges, list):
                range_lengths = [len(x) for x in ranges]
            else:
                range_lengths = [len(ranges)]
        self._range_lengths = _sanitize_range_lengths(range_lengths)

        if mcols is None:
            mcols = BiocFrame(number_of_rows=len(range_lengths))
        self._mcols = mcols

        if names is not None and not isinstance(names, ut.Names):
            names = ut.Names(names)
        self._names = names

        self._metadata = {} if metadata is None else metadata

        if validate is True:
            _validate_ranges(self._ranges, sum(self._range_lengths))
            _validate_optional_attrs(self._mcols, self._names, len(self._range_lengths))

    def _define_output(self, in_place: bool = False) -> "GenomicRangesList":
        if in_place is True:
            return self
        else:
            return self.__copy__()

    #########################
    ######>> Copying <<######
    #########################

    def __deepcopy__(self, memo=None, _nil=[]):
        """
        Returns:
            A deep copy of the current ``GenomicRangesList``.
        """
        from copy import deepcopy

        _ranges_copy = deepcopy(self._ranges)
        _rangelengths_copy = deepcopy(self._range_lengths)
        _names_copy = deepcopy(self._names)
        _mcols_copy = deepcopy(self._mcols)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            ranges=_ranges_copy,
            range_lengths=_rangelengths_copy,
            names=_names_copy,
            mcols=_mcols_copy,
            metadata=_metadata_copy,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``GenomicRangesList``.
        """
        current_class_const = type(self)
        return current_class_const(
            ranges=self._ranges,
            range_lengths=self._range_lengths,
            names=self._names,
            mcols=self._mcols,
            metadata=self._metadata,
        )

    def copy(self):
        """Alias for :py:meth:`~__copy__`."""
        return self.__copy__()

    ######################################
    ######>> length and iterators <<######
    ######################################

    def __len__(self) -> int:
        """
        Returns:
            Number of rows.
        """
        return len(self._range_lengths)

    def __iter__(self) -> GenomicRangesListIter:
        """Iterator over rows."""
        return GenomicRangesListIter(self)

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation of this ``GenomicRangesList``.
        """
        output = "GenomicRangesList(number_of_elements=" + str(len(self))
        output += ", ranges=" + repr(self._ranges)

        if self._names is not None:
            output += ", names=" + ut.print_truncated_list(self._names)

        if self._mcols is not None:
            output += ", mcols=" + repr(self._mcols)

        if len(self._metadata) > 0:
            output += ", metadata=" + ut.print_truncated_dict(self._metadata)

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this ``GenomicRangesList``.
        """
        output = (
            f"GenomicRangesList with {len(self)} range{'s' if len(self) != 1 else ''}"
        )
        output += f" and {len(self._mcols)} metadata column{'s' if len(self._mcols) != 1 else ''}\n"

        if isinstance(self._ranges, GenomicRanges) and len(self._ranges) == 0:
            output += "--- empty genomic ranges list ---"
            return output

        output += " \n"

        nr = len(self)
        added_table = False
        if nr:
            if nr <= 10:
                indices = range(nr)
                insert_ellipsis = False
            else:
                indices = [0, 1, 2, nr - 3, nr - 2, nr - 1]
                insert_ellipsis = True

            raw_floating = ut.create_floating_names(self._names, indices)
            if insert_ellipsis:
                raw_floating = raw_floating[:3] + [""] + raw_floating[3:]
            floating = raw_floating

            for idx in range(len(indices)):
                output += f"Name: {floating[idx]} \n"
                output += self._ranges[indices[idx]].__str__()
                output += "\n \n"

        footer = []
        if self._mcols is not None and self._mcols.shape[1]:
            footer.append(
                "mcols("
                + str(self._mcols.shape[1])
                + " columns): "
                + ut.print_truncated_list(
                    self._mcols.column_names,
                    sep=" ",
                    include_brackets=False,
                    transform=lambda y: y,
                )
            )
        if len(self._metadata):
            footer.append(
                "metadata("
                + str(len(self.metadata))
                + "): "
                + ut.print_truncated_list(
                    list(self.metadata.keys()),
                    sep=" ",
                    include_brackets=False,
                    transform=lambda y: y,
                )
            )
        if len(footer):
            if added_table:
                output += "\n------\n"
            output += "\n".join(footer)

        return output

    ##########################
    ######>> ranges <<########
    ##########################

    def get_ranges(self) -> Union[GenomicRanges, List[GenomicRanges]]:
        """
        Returns:
            List of genomic ranges.
        """

        return self._ranges

    def set_ranges(
        self, ranges: Union[GenomicRanges, List[GenomicRanges]], in_place: bool = False
    ) -> "GenomicRanges":
        """Set new genomic ranges.

        Args:
            ranges:
                List of genomic elements.
                All elements in this list must be
                :py:class:`genomicranges.GenomicRanges.GenomicRanges` objects.

            in_place:
                Whether to modify the ``GenomicRangesList`` object in place.

        Returns:
            A modified ``GenomicRangesList`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _validate_ranges(ranges, sum(self._range_lengths))
        output = self._define_output(in_place)
        output._ranges = ranges
        return output

    @property
    def ranges(self) -> GenomicRanges:
        """Alias for :py:meth:`~get_ranges`."""
        return self.get_ranges()

    @ranges.setter
    def ranges(self, ranges: GenomicRanges):
        """Alias for :py:meth:`~set_ranges` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'ranges' is an in-place operation, use 'set_ranges' instead",
            UserWarning,
        )
        self.set_ranges(ranges, in_place=True)

    ########################
    ######>> names <<#######
    ########################

    def get_names(self) -> ut.Names:
        """
        Returns:
            A list of names for each genomic element.
        """
        return self._names

    def set_names(
        self,
        names: Optional[Sequence[str]],
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Set new names.

        Args:
            names:
                Names for each genomic element.

            in_place:
                Whether to modify the ``GenomicRangesList`` object in place.

        Returns:
            A modified ``GenomicRangesList`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if names is not None:
            _validate_optional_attrs(None, names, len(self))

            if not isinstance(names, ut.Names):
                names = ut.Names(names)

        output = self._define_output(in_place)
        output._names = names
        return output

    @property
    def names(self) -> ut.Names:
        """Alias for :py:meth:`~get_names`."""
        return self.get_names()

    @names.setter
    def names(
        self,
        names: Optional[Sequence[str]],
    ):
        """Alias for :py:meth:`~set_names` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'names' is an in-place operation, use 'set_names' instead",
            UserWarning,
        )
        self.set_names(names, in_place=True)

    ########################
    ######>> mcols <<#######
    ########################

    def get_mcols(self) -> BiocFrame:
        """
        Returns:
            A ~py:class:`~biocframe.BiocFrame.BiocFrame` containing
            per-genomic element annotations.
        """
        return self._mcols

    def set_mcols(
        self,
        mcols: Optional[BiocFrame],
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Set new range metadata.

        Args:
            mcols:
                A ~py:class:`~biocframe.BiocFrame.BiocFrame` with
                length same as the number of genomic elements, containing
                per-genomic element annotations.

            in_place:
                Whether to modify the ``GenomicRangesList`` object in place.

        Returns:
            A modified ``GenomicRangesList`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        if mcols is None:
            mcols = BiocFrame({}, number_of_rows=len(self))

        _validate_optional_attrs(mcols, None, len(self))

        output = self._define_output(in_place)
        output._mcols = mcols
        return output

    @property
    def mcols(self) -> BiocFrame:
        """Alias for :py:meth:`~get_mcols`."""
        return self.get_mcols()

    @mcols.setter
    def mcols(
        self,
        mcols: Optional[BiocFrame],
    ):
        """Alias for :py:meth:`~set_mcols` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'mcols' is an in-place operation, use 'set_mcols' instead",
            UserWarning,
        )
        self.set_mcols(mcols, in_place=True)

    ###########################
    ######>> metadata <<#######
    ###########################

    def get_metadata(self) -> dict:
        """
        Returns:
            Dictionary of metadata for this object.
        """
        return self._metadata

    def set_metadata(self, metadata: dict, in_place: bool = False) -> "GenomicRanges":
        """Set additional metadata.

        Args:
            metadata:
                New metadata for this object.

            in_place:
                Whether to modify the ``GenomicRangesList`` object in place.

        Returns:
            A modified ``GenomicRangesList`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if not isinstance(metadata, dict):
            raise TypeError(
                f"`metadata` must be a dictionary, provided {type(metadata)}."
            )

        output = self._define_output(in_place)
        output._metadata = metadata
        return output

    @property
    def metadata(self) -> dict:
        """Alias for :py:attr:`~get_metadata`."""
        return self.get_metadata()

    @metadata.setter
    def metadata(self, metadata: dict):
        """Alias for :py:attr:`~set_metadata` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'metadata' is an in-place operation, use 'set_metadata' instead",
            UserWarning,
        )
        self.set_metadata(metadata, in_place=True)

    ###########################
    ######>> groups <<#########
    ###########################

    def groups(self, group: Union[str, int]) -> "GenomicRangesList":
        """Get a genomic element by their name or index position.

        Args:
            group:
                Name or index of the genomic element to access.

        Returns:
            The genomic element for group `i`.
        """

        if isinstance(group, str):
            group = self._names.map(group)

        if group < 0 or group > len(self):
            raise ValueError(
                "'group' must be less than the number of genomic elements."
            )

        return self[group]

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
            An integer vector where each value corresponds to the length of
            the contained GenomicRanges object.
        """
        return self._generic_accessor("__len__", func=True)

    def is_empty(self) -> bool:
        """Whether ``GRangesList`` has no elements or if all its elements are empty.

        Returns:
            True if the object has no elements.
        """
        if len(self) == 0:
            return True

        element_lengths = self.element_nrows()
        if all([True if x == 0 else False for x in element_lengths]):
            return True

        return False

    ##############################
    ######>> accessors <<#########
    ##############################

    # TODO: convert some of these properties to a factorized array

    @property
    def seqnames(self) -> Dict[str, List[str]]:
        """Get all sequence names.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list of sequence names.
        """
        return self._generic_accessor("seqnames")

    @property
    def start(self) -> Dict[str, List[int]]:
        """Get all start positions.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("start")

    @property
    def end(self) -> Dict[str, List[int]]:
        """Get all end positions.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("end")

    @property
    def width(self) -> Dict[str, List[int]]:
        """Get width of all regions across all elements.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("width")

    @property
    def strand(self) -> Dict[str, List[int]]:
        """Get strand of all regions across all elements.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("strand")

    @property
    def seq_info(self) -> Dict[str, List[int]]:
        """Get information about the underlying sequences.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("seq_info")

    @property
    def is_circular(self) -> Dict[str, List[int]]:
        """Get the circularity flag.

        Returns:
            A list with the same length as keys in the object,
            each element in the list contains another list values.
        """
        return self._generic_accessor("is_circular")

    def get_range_lengths(self) -> dict:
        """
        Returns:
            Number of ranges for each genomic element.
        """
        return self._range_lengths

    @property
    def range_lengths(self) -> dict:
        """Alias for :py:attr:`~get_range_lengths`."""
        return self.get_range_lengths()

    ###################################
    ######>> pandas interop <<#########
    ###################################

    def to_pandas(self) -> "pandas.DataFrame":
        """Coerce object to a :py:class:`pandas.DataFrame`.

        Returns:
            A :py:class:`~pandas.DataFrame` object.
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

    ############################
    ######>> slicers <<#########
    ############################

    def __getitem__(
        self, args: Union[str, int, tuple, list, slice]
    ) -> Union[GenomicRanges, "GenomicRangesList"]:
        """Subset individual genomic elements.

        Args:
            args:
                Name of the genomic element to access.

                Alternatively, if names of genomic elements are not
                available, you may provide an index position of the
                genomic element to access.

                Alternatively, ``args`` may also specify a list of
                positions to slice specified either as a
                :py:class:`~list` or :py:class:`~slice`.

                A tuple may also be specified along each dimension.
                Currently if the tuple contains more than
                one dimension, its ignored.

        Raises:
            TypeError:
                If ``args`` is not a supported slice argument.

        Returns:
            A new ``GenomicRangesList`` of the slice.
        """
        if isinstance(args, int):
            return self._ranges[args]
        elif isinstance(args, str):
            if self.names is not None:
                _idx = self.names.map(args)
                return self._ranges[_idx]
        else:
            idx, _ = ut.normalize_subscript(args, len(self), self._names)

            if isinstance(idx, list):
                if ut.is_list_of_type(idx, bool):
                    if len(idx) != len(self):
                        raise ValueError(
                            "`indices` is a boolean vector, length should match the size of the data."
                        )

                    idx = [i for i in range(len(idx)) if idx[i] is True]

                new_ranges = [self.ranges[i] for i in idx]
                new_range_lengths = [self._range_lengths[i] for i in idx]

                new_names = None
                if self.names is not None:
                    new_names = [self.names[i] for i in idx]

                new_mcols = None
                if self.mcols is not None:
                    new_mcols = self.mcols[idx, :]

                return GenomicRangesList(
                    new_ranges, new_range_lengths, new_names, new_mcols, self._metadata
                )
            elif isinstance(idx, (slice, range)):
                if isinstance(idx, range):
                    idx = slice(idx.start, idx.stop, idx.step)

                return GenomicRangesList(
                    self._ranges[idx],
                    self._range_lengths[idx],
                    self._names[idx] if self._names is not None else self._names,
                    self._mcols[idx, :],
                    self._metadata,
                )

            raise TypeError("Arguments to subset `GenomicRangesList` is not supported.")

    ##########################
    ######>> empty <<#########
    ##########################

    @classmethod
    def empty(cls, n: int):
        """Create an empty ``n``-length `GenomicRangesList` object.

        Returns:
            same type as caller, in this case a `GenomicRangesList`.
        """
        _range_lengths = [0] * n

        return cls(ranges=GenomicRanges.empty(), range_lengths=_range_lengths)


@ut.combine_sequences.register(GenomicRangesList)
def _combine_grl(*x: GenomicRangesList):
    has_names = False
    for y in x:
        if y._names is not None:
            has_names = True
            break

    all_names = None
    if has_names:
        all_names = []
        for y in x:
            if y._names is not None:
                all_names += y._names
            else:
                all_names += [""] * len(y)

    return GenomicRangesList(
        ranges=ut.combine_sequences(*[y._ranges for y in x]),
        range_lengths=ut.combine_sequences(*[y._range_lengths for y in x]),
        names=all_names,
        mcols=ut.relaxed_combine_rows(*[y._mcols for y in x]),
        metadata=x[0]._metadata,
        validate=False,
    )
