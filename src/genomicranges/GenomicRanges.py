import math
from typing import List, Literal, Optional, Sequence, Union
from warnings import warn

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from .SeqInfo import SeqInfo, merge_SeqInfo
from .utils import sanitize_strand_vector

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _guess_num_ranges(seqnames, ranges):
    if len(seqnames) != len(ranges):
        raise ValueError("Number of 'seqnames' and 'ranges' do not match!")

    return len(seqnames)


def _validate_seqnames(seqnames, seqinfo, num_ranges):
    if seqnames is None:
        raise ValueError("'seqnames' cannot be None!")

    if len(seqnames) != num_ranges:
        raise ValueError(
            "Length of 'seqnames' does not match the number of ranges.",
            f"Need to be {num_ranges}, provided {len(seqnames)}.",
        )

    if any(x is None for x in seqnames):
        raise ValueError("'seqnames' cannot contain None values.")

    if not isinstance(seqinfo, SeqInfo):
        raise TypeError("'seqinfo' is not an instance of `SeqInfo` class.")

    if not all(x < len(seqinfo) for x in seqnames):
        raise ValueError(
            "'seqnames' contains sequence name not represented in 'seqinfo'."
        )


def _validate_ranges(ranges, num_ranges):
    if ranges is None:
        raise ValueError("'ranges' cannot be None.")

    if not isinstance(ranges, IRanges):
        raise TypeError("'ranges' is not an instance of `IRanges`.")

    if len(ranges) != num_ranges:
        raise ValueError(
            "Length of 'ranges' does not match the number of seqnames.",
            f"Need to be {num_ranges}, provided {len(ranges)}.",
        )


def _validate_optional_attrs(strand, mcols, names, num_ranges):
    if strand is not None:
        if len(strand) != num_ranges:
            raise ValueError("Length of 'strand' does not match the number of ranges.")

        if any(x is None for x in strand):
            raise ValueError("'strand' cannot contain None values.")

    if mcols is not None:
        if not isinstance(mcols, BiocFrame):
            raise TypeError("'mcols' is not a `BiocFrame` object.")

        if mcols.shape[0] != num_ranges:
            raise ValueError("Length of 'mcols' does not match the number of ranges.")

    if names is not None:
        if len(names) != num_ranges:
            raise ValueError("Length of 'names' does not match the number of ranges.")

        if any(x is None for x in names):
            raise ValueError("'strand' cannot contain None values.")


class GenomicRangesIter:
    """An iterator to a :py:class:`~GenomicRanges` object.

    Args:
        obj (GenomicRanges): Source object to iterate.
    """

    def __init__(self, obj: "GenomicRanges") -> None:
        """Initialize the iterator.

        Args:
            obj (GenomicRanges): source object to iterate.
        """
        self._gr = obj
        self._current_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._current_index < len(self._gr):
            iter_row_index = (
                self._gr.names[self._current_index]
                if self._gr.names is not None
                else None
            )

            iter_slice = self._gr[self._current_index]
            self._current_index += 1
            return (iter_row_index, iter_slice)

        raise StopIteration


class GenomicRanges:
    """``GenomicRanges`` provides a container class to represent and operate over genomic regions and annotations.

    Note: The documentation for some of the methods are derived from the
    `GenomicRanges R/Bioconductor package <https://github.com/Bioconductor/GenomicRanges>`_.
    """

    def __init__(
        self,
        seqnames: Sequence[str],
        ranges: IRanges,
        strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]] = None,
        names: Optional[Sequence[str]] = None,
        mcols: Optional[BiocFrame] = None,
        seqinfo: Optional[SeqInfo] = None,
        metadata: Optional[dict] = None,
        validate: bool = True,
    ):
        """Initialize a ``GenomicRanges`` object.

        Args:
            seqnames:
                List of sequence or chromosome names.

            ranges:
                Genomic positions and widths of each position. Must have the same length as ``seqnames``.

            strand:
                Strand information for each genomic range. This should be 0 (any strand),
                1 (forward strand) or -1 (reverse strand). If None, all genomic ranges
                are assumed to be 0.

                May be provided as a list of strings representing the strand;
                "+" for forward strand, "-" for reverse strand, or "*" for any strand and will be mapped
                accordingly to 1, -1 or 0.

            names:
                Names for each genomic range. Defaults to None, which means the ranges are unnamed.

            mcols:
                A ~py:class:`~biocframe.BiocFrame.BiocFrame` with the number of rows same as
                number of genomic ranges, containing per-range annotation. Defaults to None, in which case an empty
                BiocFrame object is created.

            seqinfo:
                Sequence information. Defaults to None, in which case a
                :py:class:`~genomicranges.SeqInfo.SeqInfo` object is created with the unique set of
                chromosome names from ``seqnames``.

            metadata:
                Additional metadata. Defaults to None, and is assigned to an empty dictionary.

            validate:
                Internal use only.
        """
        if seqinfo is None:
            seqinfo = SeqInfo(seqnames=list(set(seqnames)))
        self._seqinfo = seqinfo

        self._reverse_seqindex = None
        self._seqnames = self._sanitize_seqnames(seqnames, self._seqinfo)
        self._ranges = ranges

        if strand is None:
            strand = np.repeat(0, len(self._seqnames))
        self._strand = sanitize_strand_vector(strand)

        if names is not None and not isinstance(names, ut.Names):
            names = ut.Names(names)
        self._names = names

        if mcols is None:
            mcols = BiocFrame({}, number_of_rows=len(self._seqnames))
        self._mcols = mcols

        self._metadata = metadata if metadata is not None else {}

        if validate is True:
            _num_ranges = _guess_num_ranges(self._seqnames, self._ranges)
            _validate_ranges(self._ranges, _num_ranges)
            _validate_seqnames(self._seqnames, self._seqinfo, _num_ranges)
            _validate_optional_attrs(
                self._strand, self._mcols, self._names, _num_ranges
            )

    def _build_reverse_seqindex(self, seqinfo: SeqInfo):
        self._reverse_seqindex = ut.reverse_index.build_reverse_index(seqinfo.seqnames)

    def _remove_reverse_seqindex(self):
        del self._reverse_seqindex

    def _sanitize_seqnames(self, seqnames, seqinfo):
        if self._reverse_seqindex is None:
            self._build_reverse_seqindex(seqinfo)

        if not isinstance(seqnames, np.ndarray):
            seqnames = np.array([self._reverse_seqindex[x] for x in list(seqnames)])

        return seqnames

    def _define_output(self, in_place: bool = False) -> "GenomicRanges":
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
            A deep copy of the current ``GenomicRanges``.
        """
        from copy import deepcopy

        _ranges_copy = deepcopy(self._ranges)
        _seqnames_copy = deepcopy(self._seqnames)
        _strand_copy = deepcopy(self._strand)
        _names_copy = deepcopy(self._names)
        _mcols_copy = deepcopy(self._mcols)
        _seqinfo_copy = deepcopy(self._seqinfo)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            seqnames=_seqnames_copy,
            ranges=_ranges_copy,
            strand=_strand_copy,
            names=_names_copy,
            mcols=_mcols_copy,
            seqinfo=_seqinfo_copy,
            metadata=_metadata_copy,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``GenomicRanges``.
        """
        current_class_const = type(self)
        new_instance = current_class_const(
            seqnames=self._seqnames,
            ranges=self._ranges,
            strand=self._strand,
            names=self._names,
            mcols=self._mcols,
            seqinfo=self._seqinfo,
            metadata=self._metadata,
        )

        return new_instance

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
        return len(self._ranges)

    def __iter__(self) -> GenomicRangesIter:
        """Iterator over rows."""
        return GenomicRangesIter(self)

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation of this ``GenomicRanges``.
        """
        output = "GenomicRanges(number_of_ranges=" + str(len(self))
        output += ", seqnames=" + ut.print_truncated_list(self._seqnames)
        output += ", ranges=" + repr(self._ranges)

        if self._strand is not None:
            output += ", strand=" + ut.print_truncated_list(self._strand)

        if self._names is not None:
            output += ", names=" + ut.print_truncated_list(self._names)

        if self._mcols is not None:
            output += ", mcols=" + repr(self._mcols)

        if self._seqinfo is not None:
            output += ", seqinfo" + repr(self._seqinfo)

        if len(self._metadata) > 0:
            output += ", metadata=" + ut.print_truncated_dict(self._metadata)

        output += ")"
        return output

    # TODO: needs some work!!
    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this ``GenomicRanges``.
        """
        output = f"GenomicRanges with {len(self)} range{'s' if len(self) != 1 else ''}"
        output += f" and {len(self._mcols)} metadata column{'s' if len(self._mcols) != 1 else ''}\n"

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
            floating = ["", ""] + raw_floating

            def show_slot(data, colname):
                showed = ut.show_as_cell(data, indices)
                header = [colname, "<" + ut.print_type(data) + ">"]
                showed = ut.truncate_strings(
                    showed, width=max(40, len(header[0]), len(header[1]))
                )
                if insert_ellipsis:
                    showed = showed[:3] + ["..."] + showed[3:]

                return header + showed

            columns = []
            columns.append(show_slot(self._seqnames, "seqnames"))
            columns.append(str(self._ranges))

            if self._strand is not None:
                columns.append(show_slot(self._strand, "strand"))

            # if self._mcols is not None:
            #     for col in self._mcols.colnames:
            #         columns.append(show_slot(self._mcols[col], col))

            output += ut.print_wrapped_table(columns, floating_names=floating)
            added_table = True

        footer = []
        if self._seqinfo is not None and len(self._seqinfo):
            footer.append(
                "seqinfo("
                + str(len(self._seqinfo))
                + " sequences): "
                + ut.print_truncated_list(
                    self._seqinfo.seqnames,
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
    ######>> seqnames <<######
    ##########################

    def get_seqnames(
        self, as_type: Literal["factor", "list"] = "list"
    ) -> Union[Union[np.ndarray, List[str]], np.ndarray]:
        """
        Returns:
            List of sequence names.
        """

        if as_type == "factor":
            return self._seqnames, self._seqinfo.seqnames
        elif as_type == "list":
            return [self._seqinfo.seqnames[x] for x in self._seqnames]
        else:
            raise ValueError("Argument 'as_type' must be either 'factor' or 'list'.")

    def set_seqnames(
        self, seqnames: Union[Sequence[str], np.ndarray], in_place: bool = False
    ) -> "GenomicRanges":
        """Set new sequence names.

        Args:
            seqnames:
                List of sequence or chromosome names. Optionally can be a numpy array with indices mapped
                to :py:attr:``~seqinfo``.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _validate_seqnames(seqnames, len(self))

        if not isinstance(seqnames, np.ndarray):
            seqnames = np.array(
                [self._seqinfo.seqnames.index(x) for x in list(seqnames)]
            )

        output = self._define_output(in_place)
        output._seqnames = seqnames
        return output

    @property
    def seqnames(self) -> Union[Union[np.ndarray, List[str]], np.ndarray]:
        """Alias for :py:meth:`~get_seqnames`."""
        return self.get_seqnames()

    @seqnames.setter
    def seqnames(self, seqnames: Union[Sequence[str], np.ndarray]):
        """Alias for :py:meth:`~set_seqnames` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'seqnames' is an in-place operation, use 'set_seqnames' instead",
            UserWarning,
        )
        self.set_row_names(seqnames, in_place=True)

    ########################
    ######>> ranges <<######
    ########################

    def get_ranges(self) -> IRanges:
        """
        Returns:
            An ``IRanges`` object containing the positions.
        """

        return self._ranges

    def set_ranges(self, ranges: IRanges, in_place: bool = False) -> "GenomicRanges":
        """Set new ranges.

        Args:
            ranges:
                IRanges containing the genomic positions and widths.
                Must have the same length as ``seqnames``.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        _validate_ranges(ranges, len(self))

        output = self._define_output(in_place)
        output._ranges = ranges
        return output

    @property
    def ranges(self) -> IRanges:
        """Alias for :py:meth:`~get_ranges`."""
        return self.get_seqnames()

    @ranges.setter
    def ranges(self, ranges: IRanges):
        """Alias for :py:meth:`~set_ranges` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'ranges' is an in-place operation, use 'set_ranges' instead",
            UserWarning,
        )
        self.set_ranges(ranges, in_place=True)

    ########################
    ######>> strand <<######
    ########################

    def get_strand(self) -> np.ndarray:
        """
        Returns:
            A numpy vector representing strand, 0
            for any strand, -1 for reverse strand
            and 1 for forward strand.
        """
        return self._strand

    def set_strand(
        self,
        strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]],
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Set new strand information.

        Args:
            strand:
                Strand information for each genomic range. This should be 0 (any strand),
                1 (forward strand) or -1 (reverse strand). If None, all genomic ranges
                are assumed to be 0.

                May be provided as a list of strings representing the strand;
                "+" for forward strand, "-" for reverse strand, or "*" for any strand and will be mapped
                accordingly to 1, -1 or 0.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if strand is None:
            strand = np.repeat(0, len(self))

        strand = sanitize_strand_vector(strand)
        _validate_optional_attrs(strand, None, None, len(self))
        output = self._define_output(in_place)
        output._strand = strand
        return output

    @property
    def strand(self) -> np.ndarray:
        """Alias for :py:meth:`~get_strand`."""
        return self.get_strand()

    @strand.setter
    def strand(
        self,
        strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]],
    ):
        """Alias for :py:meth:`~set_strand` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'strand' is an in-place operation, use 'set_strand' instead",
            UserWarning,
        )
        self.set_strand(strand, in_place=True)

    ########################
    ######>> names <<#######
    ########################

    def get_names(self) -> ut.Names:
        """
        Returns:
            A list of names for each genomic range.
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
                Names for each genomic range.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if names is not None:
            _validate_optional_attrs(None, None, names, len(self))

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
            A ~py:class:`~biocframe.BiocFrame.BiocFrame` containing per-range annotations.
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
                A ~py:class:`~biocframe.BiocFrame.BiocFrame` with length same as the number
                of ranges, containing per-range annotations.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        if mcols is None:
            mcols = BiocFrame({}, number_of_rows=len(self))

        _validate_optional_attrs(None, mcols, None, len(self))

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

    ##########################
    ######>> seqinfo <<#######
    ##########################

    def get_seqinfo(self) -> SeqInfo:
        """
        Returns:
            A ~py:class:`~genomicranges.SeqInfo.SeqInfo` containing sequence information.
        """
        return self._seqinfo

    def set_seqinfo(
        self,
        seqinfo: Optional[SeqInfo],
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Set new sequence information.

        Args:
            seqinfo:
                A ~py:class:`~genomicranges.SeqInfo.SeqInfo` object contaning information
                about sequences in :py:attr:`~seqnames`.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if seqinfo is None:
            seqinfo = SeqInfo(seqnames=self.get_seqnames(as_type="list"))

        _validate_seqnames(self.seqnames, seqinfo)

        output = self._define_output(in_place)
        output._seqinfo = seqinfo
        return output

    @property
    def seqinfo(self) -> np.ndarray:
        """Alias for :py:meth:`~get_seqinfo`."""
        return self.get_seqinfo()

    @seqinfo.setter
    def seqinfo(
        self,
        seqinfo: Optional[SeqInfo],
    ):
        """Alias for :py:meth:`~set_seqinfo` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'seqinfo' is an in-place operation, use 'set_seqinfo' instead",
            UserWarning,
        )
        self.set_seqinfo(seqinfo, in_place=True)

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
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
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

    ################################
    ######>> Single getters <<######
    ################################

    def get_start(self) -> np.ndarray:
        """Get all start positions.

        Returns:
            NumPy array of 32-bit signed integers containing the start
            positions for all ranges.
        """
        return self._ranges.get_start()

    @property
    def start(self) -> np.ndarray:
        """Alias for :py:attr:`~get_start`."""
        return self.get_start()

    def get_end(self) -> np.ndarray:
        """Get all end positions.

        Returns:
            NumPy array of 32-bit signed integers containing the end
            positions for all ranges.
        """
        return self._ranges.get_end()

    @property
    def end(self) -> np.ndarray:
        """Alias for :py:attr:`~get_end`."""
        return self.get_end()

    def get_width(self) -> np.ndarray:
        """Get width of each genomic range.

        Returns:
            NumPy array of 32-bit signed integers containing the width
            for all ranges.
        """
        return self._ranges.get_width()

    @property
    def width(self) -> np.ndarray:
        """Alias for :py:attr:`~get_width`."""
        return self.get_width()

    #########################
    ######>> Slicers <<######
    #########################

    def get_subset(self, subset: Union[str, int, bool, Sequence]) -> "GenomicRanges":
        """Subset ``GenomicRanges``, based on their indices or names.

        Args:
            subset:
                Indices to be extracted. This may be an integer, boolean, string,
                or any sequence thereof, as supported by
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.
                Scalars are treated as length-1 sequences.

                Strings may only be used if :py:attr:``~names`` are available (see
                :py:meth:`~get_names`). The first occurrence of each string
                in the names is used for extraction.

        Returns:
            A new ``GenomicRanges`` object with the ranges of interest.
        """
        idx, _ = ut.normalize_subscript(subset, len(self), self._names)

        current_class_const = type(self)
        return current_class_const(
            seqnames=[self._seqinfo.seqnames[x] for x in self._seqnames[idx]],
            ranges=self._ranges[idx],
            strand=self._strand[idx],
            names=self._names[idx] if self._names is not None else None,
            mcols=self._mcols[idx, :],
            seqinfo=self._seqinfo,
            metadata=self._metadata,
        )

    def __getitem__(self, subset: Union[str, int, bool, Sequence]) -> "GenomicRanges":
        """Alias to :py:attr:`~get_subset`."""
        return self.get_subset(subset)

    def set_subset(
        self,
        args: Union[Sequence, int, str, bool, slice, range],
        value: "GenomicRanges",
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Add or update positions (in-place operation).

        Args:
            args:
                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be replaced, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

            value:
                An ``GenomicRanges`` object of length equal to the number of ranges
                to be replaced, as specified by ``subset``.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        idx, _ = ut.normalize_subscript(args, len(self), self._names)

        output = self._define_output(in_place)

        output._seqnames[idx] = value._seqnames
        output._ranges[idx] = value._ranges
        output._strand[idx] = value._strand

        if value._names is not None:
            if output._names is None:
                output._names = [""] * len(output)
            for i, j in enumerate(idx):
                output._names[j] = value._names[i]
        elif output._names is not None:
            for i, j in enumerate(idx):
                output._names[j] = ""

        if value._mcols is not None:
            output._mcols[idx, :] = value._mcols

        return output

    def __setitem__(
        self,
        args: Union[Sequence, int, str, bool, slice, range],
        value: "GenomicRanges",
    ) -> "GenomicRanges":
        """Alias to :py:attr:`~set_subset`.

        This operation modifies object in-place.
        """

        warn(
            "Modifying a subset of the object is an in-place operation, use 'set_subset' instead",
            UserWarning,
        )

        return self.set_subset(args, value, in_place=True)

    ################################
    ######>> pandas interop <<######
    ################################

    def to_pandas(self) -> "pandas.DataFrame":
        """Convert this ``GenomicRanges`` object into a :py:class:`~pandas.DataFrame`.

        Returns:
            A :py:class:`~pandas.DataFrame` object.
        """
        import pandas as pd

        _rdf = self._ranges.to_pandas()
        _rdf["seqnames"] = self._seqnames
        _rdf["strand"] = self._strand

        if self._names is not None:
            _rdf.index = self._names

        if self._mcols is not None:
            if self._mcols.shape[1] > 0:
                _rdf = pd.concat([_rdf, self._mcols.to_pandas()], axis=1)

        return _rdf

    @classmethod
    def from_pandas(cls, input: "pandas.DataFrame") -> "GenomicRanges":
        """Create a ``GenomicRanges`` from a :py:class:`~pandas.DataFrame` object.

        Args:
            input:
                Input data. must contain columns 'seqnames', 'start' and 'width'.

        Returns:
            A ``GenomicRanges`` object.
        """
        from pandas import DataFrame

        if not isinstance(input, DataFrame):
            raise TypeError("`input` is not a pandas `DataFrame` object.")

        if "start" not in input.columns:
            raise ValueError("'input' must contain column 'start'.")
        start = input["start"].tolist()

        if "width" not in input.columns:
            raise ValueError("'input' must contain column 'width'.")
        width = input["width"].tolist()

        if "seqnames" not in input.columns:
            raise ValueError("'input' must contain column 'seqnames'.")
        seqnames = input["seqnames"].tolist()

        ranges = IRanges(start, width)

        # mcols
        mcols_df = input.drop(columns=["start", "width", "seqnames"])

        mcols = None
        if (not mcols_df.empty) or len(mcols_df.columns) > 0:
            mcols = BiocFrame.from_pandas(mcols_df)

        names = None
        if input.index is not None:
            names = [str(i) for i in input.index.to_list()]

        return cls(ranges=ranges, seqnames=seqnames, names=names, mcols=mcols)

    #####################################
    ######>> intra-range methods <<######
    #####################################

    def flank(
        self,
        width: int,
        start: bool = True,
        both: bool = False,
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Compute flanking ranges for each range. The logic for this comes from the `R/GenomicRanges` & `IRanges`
        packages.

        If ``start`` is ``True`` for a given range, the flanking occurs at the `start`,
        otherwise the `end`.
        The `widths` of the flanks are given by the ``width`` parameter.

        ``width`` can be negative, in which case the flanking region is
        reversed so that it represents a prefix or suffix of the range.

        Usage:

            `gr.flank(3, True)`, where "x" indicates a range in ``gr`` and "-" indicates the
            resulting flanking region:
                ---xxxxxxx
            If ``start`` were ``False``, the range in ``gr`` becomes
                xxxxxxx---
            For negative width, i.e. `gr.flank(x, -3, FALSE)`, where "*" indicates the overlap
            between "x" and the result:
                xxxx***
            If ``both`` is ``True``, then, for all ranges in "x", the flanking regions are
            extended into (or out of, if ``width`` is negative) the range, so that the result
            straddles the given endpoint and has twice the width given by width.

            This is illustrated below for `gr.flank(3, both=TRUE)`:
                ---***xxxx

        Args:
            width:
                Width to flank by. May be negative.

            start:
                Whether to only flank starts. Defaults to True.

            both:
                Whether to flank both starts and ends. Defaults to False.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the flanked regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """

        all_starts = self.start
        all_ends = self.end
        all_strands = self.strand

        # figure out which position to pin, start or end?
        start_flags = np.repeat(start, len(all_strands))
        if not ignore_strand:
            start_flags = [
                start != (all_strands[i] == -1) for i in range(len(all_strands))
            ]

        new_starts = []
        new_widths = []
        # if both is true, then depending on start_flag, we extend it out
        # I couldn't understand the scenario's with witdh <=0,
        # so refer to the R implementation here
        for idx in range(len(start_flags)):
            sf = start_flags[idx]
            tstart = 0
            if both is True:
                tstart = (
                    all_starts[idx] - abs(width) if sf else all_ends[idx] - abs(width)
                )
            else:
                if width >= 0:
                    tstart = all_starts[idx] - abs(width) if sf else all_ends[idx]
                else:
                    tstart = all_starts[idx] if sf else all_ends[idx] + width

            new_starts.append(tstart)
            new_widths.append((width * (2 if both else 1)))

        output = self._define_output(in_place)
        output._ranges = IRanges(new_starts, new_widths)
        return output

    def resize(
        self,
        width: int,
        fix: Literal["start", "end", "center"] = "start",
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Resize ranges to the specified ``width`` where either the ``start``, ``end``,
        or ``center`` is used as an anchor.

        Args:
            width:
                Width to resize, cannot be negative!

            fix:
                Fix positions by "start", "end", or "center".
                Defaults to "start".

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Raises:
            ValueError:
                If parameter ``fix`` is neither `start`, `end`, nor `center`.
                If ``width`` is negative.

        Returns:
            A modified ``GenomicRanges`` object with the resized regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """

        if width < 0:
            raise ValueError("`width` cannot be negative!")

        if fix not in ["start", "end", "center"]:
            raise ValueError(
                f"`fix` must be either 'start', 'end' or 'center', provided {fix}"
            )

        new_starts = []
        new_widths = []

        for idx, row in self:
            ts = None
            te = None

            _strand = row.strand[0]
            _start = row.start[0]
            _end = row.end[0]

            if ignore_strand is True or _strand != -1:
                if fix == "start":
                    ts = _start
                    te = _start + width - 1
                elif fix == "center":
                    tmid = math.ceil((_start + _end) / 2)
                    twidthby2 = (
                        math.floor(width / 2) if _strand == 1 else math.ceil(width / 2)
                    )
                    ts = tmid - twidthby2
                    te = ts + width - 1
                else:
                    te = _end
                    ts = _end - width   
            elif _strand == -1:
                if fix == "end":
                    ts = _start
                    te = _start + width - 1
                elif fix == "center":
                    tmid = math.ceil((_start + _end) / 2)
                    twidthby2 = math.ceil(width / 2)
                    ts = tmid - twidthby2
                    te = ts + width - 1
                else:
                    te = _end
                    ts = _end - width

            new_starts.append(ts)
            new_widths.append(te - ts)

        output = self._define_output(in_place)
        output._ranges = IRanges(new_starts, new_widths)
        return output


@ut.combine_sequences.register
def _combine_GenomicRanges(*x: GenomicRanges) -> GenomicRanges:
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

    return GenomicRanges(
        ranges=ut.combine_sequences(*[y._ranges for y in x]),
        seqnames=ut.combine_sequences(*[y._seqnames for y in x]),
        strand=ut.combine_sequences(*[y._strand for y in x]),
        names=all_names,
        mcols=ut.combine_rows(*[y._mcols for y in x]),
        seqinfo=merge_SeqInfo([y._seqinfo for y in x]),
        metadata=x[0]._metadata,
        validate=False,
    )
