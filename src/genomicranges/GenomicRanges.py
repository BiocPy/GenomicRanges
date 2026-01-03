from __future__ import annotations

from collections import defaultdict
from multiprocessing import Pool, cpu_count
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple, Union
from warnings import warn

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from .sequence_info import SeqInfo, merge_SeqInfo
from .utils import (
    STRAND_MAP,
    compute_up_down,
    group_by_indices,
    sanitize_strand_vector,
    wrapper_follow_precede,
)

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

_granges_delim = "__"


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

    if np.isnan(seqnames).any():
        raise ValueError("'seqnames' cannot contain None values.")

    if not isinstance(seqinfo, SeqInfo):
        raise TypeError("'seqinfo' is not an instance of `SeqInfo` class.")

    _l = len(seqinfo)
    if (seqnames > _l).any():
        raise ValueError("'seqnames' contains sequence name not represented in 'seqinfo'.")


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

        if np.isnan(strand).any():
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
            raise ValueError("'names' cannot contain None values.")


class GenomicRangesIter:
    """An iterator to a :py:class:`~GenomicRanges` object."""

    def __init__(self, obj: GenomicRanges) -> None:
        """Initialize the iterator.

        Args:
            obj:
                Source object to iterate.
        """
        self._gr = obj
        self._current_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._current_index < len(self._gr):
            iter_row_index = self._gr.names[self._current_index] if self._gr.names is not None else None

            iter_slice = self._gr[self._current_index]
            self._current_index += 1
            return (iter_row_index, iter_slice)

        raise StopIteration


class GenomicRanges(ut.BiocObject):
    """``GenomicRanges`` provides a container class to represent and operate over genomic regions and annotations.

    Note: The documentation for some of the methods are derived from the
    `GenomicRanges R/Bioconductor package <https://github.com/Bioconductor/GenomicRanges>`_.
    """

    def __init__(
        self,
        seqnames: Sequence[str],
        ranges: IRanges,
        strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]] = None,
        names: Optional[Union[ut.Names, Sequence[str]]] = None,
        mcols: Optional[BiocFrame] = None,
        seqinfo: Optional[SeqInfo] = None,
        metadata: Optional[Union[Dict[str, Any], ut.NamedList]] = None,
        _validate: bool = True,
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
                are assumed to be 0 (any strand).

                May be provided as a list of strings representing the strand;
                "+" for forward strand, "-" for reverse strand, or "*" for any strand and will be mapped
                to 1, -1 or 0 respectively.

            names:
                Names for each genomic range.
                Defaults to None, which means the ranges are unnamed.

            mcols:
                A ~py:class:`~biocframe.BiocFrame.BiocFrame` with the number of rows same as
                number of genomic ranges, containing per-range annotation.
                Defaults to None, in which case an empty BiocFrame object is created.

            seqinfo:
                Sequence information. Defaults to None, in which case a
                :py:class:`~genomicranges.SeqInfo.SeqInfo` object is created with the unique set of
                chromosome names from ``seqnames``.

            metadata:
                Additional metadata.
                Defaults to None, and is assigned to an empty dictionary.

            validate:
                Internal use only.
        """

        super().__init__(metadata=metadata, _validate=_validate)

        if seqinfo is None:
            seqinfo = SeqInfo(seqnames=sorted(list(set(seqnames))))
        self._seqinfo = seqinfo

        self._reverse_seqindex = None
        self._seqnames = self._sanitize_seqnames(seqnames, self._seqinfo)
        self._ranges = ranges

        if strand is None:
            strand = np.repeat(0, len(self._seqnames))
            self._strand = np.asarray(strand, dtype=np.int8)
        else:
            self._strand = sanitize_strand_vector(strand)

        if names is not None and not isinstance(names, ut.Names):
            names = ut.Names(names)
        self._names = names

        if mcols is None:
            mcols = BiocFrame({}, number_of_rows=len(self._seqnames))
        self._mcols = mcols

        if _validate is True:
            _num_ranges = _guess_num_ranges(self._seqnames, self._ranges)
            _validate_ranges(self._ranges, _num_ranges)
            _validate_seqnames(self._seqnames, self._seqinfo, _num_ranges)
            _validate_optional_attrs(self._strand, self._mcols, self._names, _num_ranges)

    def _build_reverse_seqindex(self, seqinfo: SeqInfo):
        self._reverse_seqindex = ut.reverse_index.build_reverse_index(seqinfo.get_seqnames())

    def _remove_reverse_seqindex(self):
        del self._reverse_seqindex

    def _sanitize_seqnames(self, seqnames, seqinfo):
        if self._reverse_seqindex is None:
            self._build_reverse_seqindex(seqinfo)

        if not isinstance(seqnames, np.ndarray):
            seqnames = np.asarray([self._reverse_seqindex[x] for x in seqnames])

            if len(seqnames) == 0:
                seqnames = seqnames.astype(np.uint8)
            else:
                num_uniq = np.max(seqnames)
                _types = [np.uint8, np.uint16, np.uint32, np.uint64]
                for _dtype in _types:
                    if num_uniq < np.iinfo(_dtype).max:
                        seqnames = seqnames.astype(_dtype)
                        break

        return seqnames

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
            _validate=False,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``GenomicRanges``.
        """
        current_class_const = type(self)
        return current_class_const(
            seqnames=self._seqnames,
            ranges=self._ranges,
            strand=self._strand,
            names=self._names,
            mcols=self._mcols,
            seqinfo=self._seqinfo,
            metadata=self._metadata,
            _validate=False,
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

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this ``GenomicRanges``.
        """
        output = f"GenomicRanges with {len(self)} range{'s' if len(self) != 1 else ''}"
        output += f" and {len(self._mcols.get_column_names())} metadata column{'s' if len(self._mcols.get_column_names()) != 1 else ''}\n"

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

            columns = []

            header = ["seqnames", "<str>"]
            _seqnames = []
            for x in self._seqnames[indices]:
                _seqnames.append(self._seqinfo.get_seqnames()[x])

            showed = _seqnames
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            header = ["ranges", "<IRanges>"]
            showed = [f"{x._start[0]} - {x.end[0]}" for _, x in self._ranges[indices]]
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            header = ["strand", f"<{ut.print_type(self._strand)}>"]
            _strand = []
            for x in self._strand[indices]:
                if x == 1:
                    _strand.append("+")
                elif x == -1:
                    _strand.append("-")
                else:
                    _strand.append("*")

            showed = _strand
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            if self._mcols.shape[1] > 0:
                spacer = ["|"] * (len(indices) + insert_ellipsis)
                columns.append(["", ""] + spacer)

                for col in self._mcols.get_column_names():
                    data = self._mcols.column(col)
                    showed = ut.show_as_cell(data, indices)
                    header = [col, "<" + ut.print_type(data) + ">"]
                    showed = ut.truncate_strings(showed, width=max(40, len(header[0]), len(header[1])))
                    if insert_ellipsis:
                        showed = showed[:3] + ["..."] + showed[3:]
                    columns.append(header + showed)

            output += ut.print_wrapped_table(columns, floating_names=floating)
            added_table = True

        footer = []
        if self._seqinfo is not None and len(self._seqinfo):
            footer.append(
                "seqinfo("
                + str(len(self._seqinfo))
                + " sequences): "
                + ut.print_truncated_list(
                    self._seqinfo.get_seqnames(),
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

    def get_seqnames(self, as_type: Literal["factor", "list"] = "list") -> Union[ut.Factor, List[str]]:
        """Access sequence names.

        Args:
            as_type:
                Access seqnames as factor tuple, in which case, levels and codes
                are returned.

                If ``list``, then codes are mapped to levels and returned.

        Returns:
            A :py:class:`biocutils.Factor` if `as_type="factor"`.
            Otherwise a list of sequence names.
        """

        if as_type == "factor":
            return ut.Factor(codes=self._seqnames, levels=self._seqinfo.get_seqnames())
        elif as_type == "list":
            return [self._seqinfo.get_seqnames()[x] for x in self._seqnames]
        else:
            raise ValueError("Argument 'as_type' must be 'factor' or 'list'.")

    def set_seqnames(self, seqnames: Union[Sequence[str], np.ndarray], in_place: bool = False) -> GenomicRanges:
        """Set new sequence names.

        Args:
            seqnames:
                List of sequence or chromosome names.
                Optionally can be a numpy array with indices mapped
                to :py:attr:``~seqinfo``.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _validate_seqnames(seqnames, len(self))

        if not isinstance(seqnames, np.ndarray):
            seqnames = np.asarray([self._seqinfo.get_seqnames().index(x) for x in list(seqnames)])

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
        self.set_seqnames(seqnames, in_place=True)

    ########################
    ######>> ranges <<######
    ########################

    def get_ranges(self) -> IRanges:
        """
        Returns:
            An ``IRanges`` object containing the positions.
        """

        return self._ranges

    def set_ranges(self, ranges: IRanges, in_place: bool = False) -> GenomicRanges:
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
        return self.get_ranges()

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

    def get_strand(
        self, as_type: Literal["numpy", "factor", "list"] = "numpy"
    ) -> Union[Tuple[np.ndarray, dict], List[str]]:
        """Access strand information.

        Args:
            as_type:
                Access seqnames as factor codes, in which case, a numpy
                 vector is retuned.

                If ``factor``, a tuple with codes as the strand vector
                and levels a dictionary containing the mapping.

                If ``list``, then codes are mapped to levels and returned.

        Returns:
            A numpy vector representing strand, 0
            for any strand, -1 for reverse strand
            and 1 for forward strand.

            A tuple of codes and levels.

            A list of "+", "-", or "*" for each range.
        """
        if as_type == "numpy":
            return self._strand
        elif as_type == "factor":
            return self._strand, {"-1": "-", "0": "*", "1": "+"}
        elif as_type == "list":
            _strand = []
            for x in self._strand:
                if x == 1:
                    _strand.append("+")
                elif x == -1:
                    _strand.append("-")
                else:
                    _strand.append("*")
            return _strand
        else:
            raise ValueError("Argument 'as_type' must be 'factor' or 'list'.")

    def set_strand(
        self, strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]], in_place: bool = False
    ) -> GenomicRanges:
        """Set new strand information.

        Args:
            strand:
                Strand information for each genomic range. This should be 0 (any strand),
                1 (forward strand) or -1 (reverse strand).

                Alternatively, may provide a list of strings representing the strand;
                "+" for forward strand, "-" for reverse strand, or "*" for any strand
                and will be automatically mapped to 1, -1 or 0 respectively.

                May be set to None; in which case all genomic ranges are assumed to be 0.

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
    def strand(self, strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]]):
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

    def set_names(self, names: Optional[Union[ut.Names, Sequence[str]]], in_place: bool = False) -> GenomicRanges:
        """Set new names.

        Args:
            names:
                Names for each genomic range.
                May be `None` to remove names.

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
    def names(self, names: Optional[Union[ut.Names, Sequence[str]]]):
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

    def set_mcols(self, mcols: Optional[BiocFrame], in_place: bool = False) -> GenomicRanges:
        """Set new range metadata.

        Args:
            mcols:
                A ~py:class:`~biocframe.BiocFrame.BiocFrame` with length same as the number
                of ranges, containing per-range annotations.

                May be None to remove range metadata.

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
    def mcols(self, mcols: Optional[BiocFrame]):
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

    def set_seqinfo(self, seqinfo: Optional[SeqInfo], in_place: bool = False) -> GenomicRanges:
        """Set new sequence information.

        Args:
            seqinfo:
                A ~py:class:`~genomicranges.SeqInfo.SeqInfo` object containing information
                about sequences in :py:attr:`~seqnames`.

                May be None to remove sequence information. This would then generate a
                new sequence information based on the current sequence names.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        seqnames = self.get_seqnames(as_type="list")

        if seqinfo is None:
            seqinfo = SeqInfo(seqnames=seqnames)

        self._reverse_seqindex = None
        new_seqnames = self._sanitize_seqnames(seqnames, seqinfo)

        _validate_seqnames(new_seqnames, seqinfo, len(self))

        output = self._define_output(in_place)
        output._seqinfo = seqinfo
        output._seqnames = new_seqnames
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

    def get_seqlengths(self) -> np.ndarray:
        """Get sequence lengths for each genomic range.

        Returns:
            An ndarray containing the sequence lengths.
        """
        seqlenths = self._seqinfo.get_seqlengths()
        return np.asarray([seqlenths[x] for x in self._seqnames])

    #########################
    ######>> Slicers <<######
    #########################

    def get_subset(self, subset: Union[str, int, bool, Sequence]) -> GenomicRanges:
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
        if len(self) == 0:
            return GenomicRanges.empty()

        idx, _ = ut.normalize_subscript(subset, len(self), self._names)

        current_class_const = type(self)
        return current_class_const(
            seqnames=ut.subset_sequence(self._seqnames, idx),
            ranges=ut.subset_sequence(self._ranges, idx),
            strand=ut.subset_sequence(self._strand, idx),
            names=ut.subset_sequence(self._names, idx) if self._names is not None else None,
            mcols=ut.subset(self._mcols, idx),
            seqinfo=self._seqinfo,
            metadata=self._metadata,
        )

    def __getitem__(self, subset: Union[str, int, bool, Sequence]) -> GenomicRanges:
        """Alias to :py:attr:`~get_subset`."""
        return self.get_subset(subset)

    def set_subset(
        self,
        args: Union[Sequence, int, str, bool, slice, range],
        value: GenomicRanges,
        in_place: bool = False,
    ) -> GenomicRanges:
        """Udate positions.

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

        if not isinstance(value, GenomicRanges):
            raise TypeError("'value' to assign must be a `GenomicRanges` object.")

        idx, _ = ut.normalize_subscript(args, len(self), self._names)

        output = self._define_output(in_place)

        output._seqnames[idx] = value._seqnames
        output._ranges[idx] = value._ranges
        output._strand[idx] = value._strand

        if self._names is None and value._names is None:
            self._names = None
        else:
            _names = list(self._names) if self._names is not None else None
            if value._names is not None:
                if _names is None:
                    _names = [""] * len(output)
                for i, j in enumerate(idx):
                    _names[j] = value._names[i]
            elif _names is not None:
                for i, j in enumerate(idx):
                    _names[j] = ""
            self._names = ut.Names(_names)

        if value._mcols is not None:
            output._mcols[idx, :] = value._mcols

        if in_place is True:
            output._ranges.delete_nclist_index()

        return output

    def __setitem__(
        self,
        args: Union[Sequence, int, str, bool, slice, range],
        value: GenomicRanges,
    ) -> GenomicRanges:
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

    def to_pandas(self):
        """Convert this ``GenomicRanges`` to a :py:class:`~pandas.DataFrame` object.

        Returns:
            A :py:class:`~pandas.DataFrame` object.
        """
        import pandas as pd

        _rdf = self._ranges.to_pandas()
        _rdf["seqnames"] = self.get_seqnames()
        _rdf["strand"] = self.get_strand(as_type="list")

        if self._names is not None:
            _rdf.index = self._names

        if self._mcols is not None:
            if self._mcols.shape[1] > 0:
                _rdf = pd.concat([_rdf, self._mcols.to_pandas()], axis=1)

        return _rdf

    @classmethod
    def from_pandas(cls, input) -> GenomicRanges:
        """Create an ``GenomicRanges`` object from a :py:class:`~pandas.DataFrame`.

        Args:
            input:
                Input data. Must contain columns 'seqnames', 'starts' and 'widths' or "ends".
                'ends' are expected to be inclusive.

        Returns:
            A ``GenomicRanges`` object.
        """
        from pandas import DataFrame

        if not isinstance(input, DataFrame):
            raise TypeError("`input` is not a pandas `DataFrame` object.")

        if "starts" not in input.columns:
            raise ValueError("'input' must contain column 'starts'.")
        start = input["starts"].tolist()

        if "widths" not in input.columns and "ends" not in input.columns:
            raise ValueError("'input' must contain either 'widths' or 'ends' columns.")

        drops = []
        if "widths" in input.columns:
            drops.append("widths")
            width = input["widths"].tolist()
        else:
            drops.append("ends")
            width = input["ends"] - input["starts"] + 1

        if "seqnames" not in input.columns:
            raise ValueError("'input' must contain column 'seqnames'.")
        seqnames = input["seqnames"].tolist()

        ranges = IRanges(start, width)

        strand = None
        if "strand" in input.columns:
            strand = input["strand"].tolist()
            drops.append("strand")

        # mcols
        drops.extend(["starts", "seqnames"])
        mcols_df = input.drop(columns=drops)

        mcols = None
        if (not mcols_df.empty) or len(mcols_df.columns) > 0:
            mcols = BiocFrame.from_pandas(mcols_df)

        names = None
        if input.index is not None:
            names = [str(i) for i in input.index.to_list()]

        return cls(ranges=ranges, seqnames=seqnames, strand=strand, names=names, mcols=mcols)

    ################################
    ######>> polars interop <<######
    ################################

    def to_polars(self):
        """Convert this ``GenomicRanges`` to a :py:class:`~polars.DataFrame` object.

        Returns:
            A :py:class:`~polars.DataFrame` object.
        """
        import polars as pl

        _rdf = self._ranges.to_polars()
        _rdf = _rdf.with_columns(seqnames=self.get_seqnames(), strand=self.get_strand(as_type="list"))

        if self._names is not None:
            _rdf = _rdf.with_columns(rownames=self._names)

        if self._mcols is not None:
            if self._mcols.shape[1] > 0:
                _rdf = pl.concat([_rdf, self._mcols.to_polars()], how="horizontal")

        return _rdf

    @classmethod
    def from_polars(cls, input) -> GenomicRanges:
        """Create an ``GenomicRanges`` object from a :py:class:`~polars.DataFrame`.

        Args:
            input:
                Input polars DataFrame. Must contain columns 'seqnames', 'starts' and 'widths' or "ends".
                'ends' are expected to be inclusive.

        Returns:
            A ``GenomicRanges`` object.
        """
        from polars import DataFrame

        if not isinstance(input, DataFrame):
            raise TypeError("`input` is not a polars `DataFrame` object.")

        if "starts" not in input.columns:
            raise ValueError("'input' must contain column 'starts'.")
        start = input["starts"].to_list()

        if "widths" not in input.columns and "ends" not in input.columns:
            raise ValueError("'input' must contain either 'widths' or 'ends' columns.")

        drops = []
        if "widths" in input.columns:
            drops.append("widths")
            width = input["widths"].to_list()
        else:
            drops.append("ends")
            width = input["ends"] - input["starts"] + 1

        if "seqnames" not in input.columns:
            raise ValueError("'input' must contain column 'seqnames'.")
        seqnames = input["seqnames"].to_list()

        ranges = IRanges(start, width)

        strand = None
        if "strand" in input.columns:
            strand = input["strand"].to_list()
            drops.append("strand")

        # mcols
        drops.extend(["starts", "seqnames"])
        mcols_df = input.drop(drops)

        mcols = None
        if (not mcols_df.is_empty()) or len(mcols_df.columns) > 0:
            mcols = BiocFrame.from_polars(mcols_df)

        names = None

        return cls(ranges=ranges, seqnames=seqnames, strand=strand, names=names, mcols=mcols)

    #####################################
    ######>> intra-range methods <<######
    #####################################

    def flank(
        self,
        width: int,
        start: Union[bool, np.ndarray, List[bool]] = True,
        both: bool = False,
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> GenomicRanges:
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
                Whether to only flank starts.
                Defaults to True.

                Alternatively, you may provide a list of start values,
                whose length is the same as the number of ranges.

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

        if not isinstance(ignore_strand, bool):
            raise TypeError("'ignore_strand' must be True or False.")

        output = self._define_output(in_place)

        if isinstance(start, bool):
            start_arr = np.full(len(output), start, dtype=bool)
        else:
            start_arr = np.asarray(start, dtype=bool)
            if len(start_arr) != len(output):
                start_arr = np.resize(start_arr, len(output))  # may be throw an error?

        if not ignore_strand:
            start_arr = np.asarray(start_arr != (output.get_strand() == -1))

        new_ranges = output._ranges.flank(width=width, start=start_arr, both=both, in_place=False)

        output._ranges = new_ranges
        return output

    def resize(
        self,
        width: Union[int, List[int], np.ndarray],
        fix: Union[Literal["start", "end", "center"], List[Literal["start", "end", "center"]]] = "start",
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> GenomicRanges:
        """Resize ranges to the specified ``width`` where either the ``start``, ``end``, or ``center`` is used as an
        anchor.

        Args:
            width:
                Width to resize, cannot be negative!

            fix:
                Fix positions by "start", "end", or "center".
                Defaults to "start".

                Alternatively, may provide a list of fix positions
                for each genomic range.

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
        _REV_FIX = {"start": "end", "end": "start", "center": "center"}

        if not isinstance(ignore_strand, bool):
            raise TypeError("'ignore_strand' must be True or False.")

        output = self._define_output(in_place)

        if isinstance(fix, str):
            fix_arr = [fix] * len(output)
        else:
            fix_arr = fix

        if len(output) > 0 and (len(fix_arr) > len(output) or len(output) % len(fix_arr) != 0):
            raise ValueError("Number of ranges is not a multiple of 'fix' length")

        if ignore_strand:
            fix_arr = [fix_arr[i % len(fix_arr)] for i in range(len(output))]
        else:
            if len(output) == 0:
                fix_arr = []
            else:
                fix_arr = [_REV_FIX[f] if strand == -1 else f for f, strand in zip(fix_arr, output.get_strand())]

        output._ranges = self._ranges.resize(width=width, fix=fix_arr)
        return output

    def shift(self, shift: Union[int, List[int], np.ndarray] = 0, in_place: bool = False) -> GenomicRanges:
        """Shift all intervals.

        Args:
            shift:
                Shift interval. If shift is 0, the current
                object is returned. Defaults to 0.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the shifted regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        output = self._define_output(in_place)
        output._ranges = self._ranges.shift(shift=shift)
        return output

    def promoters(self, upstream: int = 2000, downstream: int = 200, in_place: bool = False) -> GenomicRanges:
        """Extend ranges to promoter regions.

        Generates promoter ranges relative to the transcription start site (TSS).

        Args:
            upstream:
                Number of positions to extend in the 5' direction.
                Defaults to 2000.

            downstream:
                Number of positions to extend in the 3' direction.
                Defaults to 200.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the promoter regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        output = self._define_output(in_place)

        new_starts, new_widths = compute_up_down(
            output.get_start(), output.get_end(), output.get_strand(), upstream, downstream, site="TSS"
        )

        output._ranges = IRanges(start=new_starts, width=new_widths)
        return output

    def terminators(self, upstream: int = 2000, downstream: int = 200, in_place: bool = False) -> GenomicRanges:
        """Extend ranges to termiantor regions.

        Generates terminator ranges relative to the transcription end site (TES).

        Args:
            upstream:
                Number of positions to extend in the 5' direction.
                Defaults to 2000.

            downstream:
                Number of positions to extend in the 3' direction.
                Defaults to 200.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the promoter regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        output = self._define_output(in_place)

        new_starts, new_widths = compute_up_down(
            output.get_start(), output.get_end(), output.get_strand(), upstream, downstream, site="TES"
        )

        output._ranges = IRanges(start=new_starts, width=new_widths)
        return output

    def restrict(
        self,
        start: Optional[Union[int, Dict[str, int], np.ndarray]] = None,
        end: Optional[Union[int, Dict[str, int], np.ndarray]] = None,
        keep_all_ranges: bool = False,
    ) -> GenomicRanges:
        """Restrict ranges to a given start and end positions.

        Args:
            start:
                Start position. Defaults to None.

                Alternatively may provide a dictionary
                mapping sequence names to starts, or array of starts.

            end:
                End position. Defaults to None.

                Alternatively may provide a dictionary
                mapping sequence names to starts, or array of starts.

            keep_all_ranges:
                Whether to keep intervals that do not overlap with start and end.
                Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with the restricted regions.
        """

        start_is_dict = isinstance(start, dict)
        end_is_dict = isinstance(end, dict)

        if not (start_is_dict or end_is_dict):
            # directly restrict ranges
            new_ranges = self._ranges.restrict(start=start, end=end, keep_all_ranges=keep_all_ranges)

            new_seqnames = self._seqnames
            new_strand = self.strand
            new_names = self._names
            new_mcols = self._mcols

            # Get indices of kept ranges for filtering seqnames and strand
            if not keep_all_ranges:
                keep_idx = new_ranges.get_mcols().get_column("revmap")

                new_seqnames = ut.subset_sequence(self._seqnames, keep_idx)
                new_strand = ut.subset_sequence(self._strand, keep_idx)
                new_names = ut.subset_sequence(self._names, keep_idx) if self._names is not None else None
                new_mcols = ut.subset(self._mcols, keep_idx)

            return GenomicRanges(
                seqnames=new_seqnames,
                ranges=new_ranges,
                strand=new_strand,
                names=new_names,
                mcols=new_mcols,
                seqinfo=self._seqinfo,
                metadata=self._metadata,
            )

        if start is None:
            start = {}
        elif not start_is_dict:
            start = {seq: start for seq in np.unique(self._seqinfo._seqnames)}

        if end is None:
            end = {}
        elif not end_is_dict:
            end = {seq: end for seq in np.unique(self._seqinfo._seqnames)}

        split_indices = defaultdict(list)
        seqnames = self._seqnames
        for i, seq in enumerate(seqnames):
            split_indices[seq].append(i)

        all_indices = []
        new_granges_list = []

        for seq in np.unique(seqnames):
            idx = split_indices[seq]
            seq_ranges = ut.subset_sequence(self._ranges, idx)
            seq_seqnames = ut.subset_sequence(seqnames, idx)
            seq_strand = ut.subset_sequence(self._strand, idx)
            seq_names = ut.subset_sequence(self._names, idx) if self._names is not None else None
            seq_mcols = ut.subset(self._mcols, idx)

            seq_start = start.get(self._seqinfo._seqnames[seq], None)
            seq_end = end.get(self._seqinfo._seqnames[seq], None)

            new_seq_ranges = seq_ranges.restrict(start=seq_start, end=seq_end, keep_all_ranges=keep_all_ranges)

            # Only keep metadata for non-dropped ranges; IRanges restrict returns a revmap
            if not keep_all_ranges:
                keep_idx = new_seq_ranges.get_mcols().get_column("revmap")
                new_seqnames = ut.subset_sequence(self._seqnames, keep_idx)
                new_strand = ut.subset_sequence(self._strand, keep_idx)
                new_names = ut.subset_sequence(self._names, keep_idx) if self._names is not None else None
                new_mcols = ut.subset(self._mcols, keep_idx)
                all_indices.extend([idx[i] for i, k in enumerate(keep_idx)])
            else:
                new_seqnames = seq_seqnames
                new_strand = seq_strand
                new_names = seq_names
                new_mcols = seq_mcols
                all_indices.extend(idx)

            new_granges_list.append(
                GenomicRanges(
                    seqnames=new_seqnames,
                    ranges=new_seq_ranges,
                    strand=new_strand,
                    names=new_names,
                    mcols=new_mcols,
                    seqinfo=self._seqinfo,
                    metadata=self._metadata,
                )
            )

        combined_granges = _combine_GenomicRanges(*new_granges_list)
        order = np.argsort(all_indices)
        ordered_granges = combined_granges[order]
        return ordered_granges

    def get_out_of_bound_index(self) -> np.ndarray:
        """Find indices of genomic ranges that are out of bounds.

        Returns:
            A numpy array of integer indices where ranges are out of bounds.
        """
        if len(self) == 0:
            return np.array([])

        # incase it contains NA's
        seqlevel_is_circ = [val is True for val in self._seqinfo._is_circular]
        seqlength_is_na = [slength is None for slength in self._seqinfo._seqlengths]
        seqlevel_has_bounds = [not (circ or is_na) for circ, is_na in zip(seqlevel_is_circ, seqlength_is_na)]

        out_of_bounds = []
        starts = self.get_start()
        ends = self.get_end()
        for i, seq_id in enumerate(self._seqnames):
            if seqlevel_has_bounds[seq_id] and (starts[i] < 1 or ends[i] > self._seqinfo._seqlengths[seq_id]):
                out_of_bounds.append(i)

        return np.asarray(out_of_bounds, dtype=np.int32)

    def trim(self, in_place: bool = False) -> GenomicRanges:
        """Trim sequences outside of bounds for non-circular chromosomes.

        Args:
            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A new ``GenomicRanges`` object with the trimmed regions.
        """

        if self.seqinfo is None:
            raise ValueError("Cannot trim ranges. `seqinfo` is not available.")

        seqlengths = self.seqinfo._seqlengths
        is_circular = self.seqinfo._is_circular

        if seqlengths is None:
            raise ValueError("Cannot trim ranges. `seqlengths` is not available.")

        if is_circular is None:
            warn("considering all sequences as non-circular...")

        out_of_bounds = self.get_out_of_bound_index()

        output = self._define_output(in_place=in_place)
        if len(out_of_bounds) == 0:
            return output.__copy__()

        filtered_seqs = self._seqnames[out_of_bounds]
        new_ends = self.get_seqlengths()[filtered_seqs]
        output._ranges[out_of_bounds] = output._ranges[out_of_bounds].restrict(
            start=1, end=new_ends, keep_all_ranges=True
        )
        return output

    def narrow(
        self,
        start: Optional[Union[int, List[int], np.ndarray]] = None,
        width: Optional[Union[int, List[int], np.ndarray]] = None,
        end: Optional[Union[int, List[int], np.ndarray]] = None,
        in_place: bool = False,
    ) -> GenomicRanges:
        """Narrow genomic positions by provided ``start``, ``width`` and ``end`` parameters.

        Important: these parameters are relative shift in positions for each range.

        Args:
            start:
                Relative start position. Defaults to None.

            width:
                Relative end position. Defaults to None.

            end:
                Relative width of the interval. Defaults to None.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the trimmed regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        if start is not None and end is not None and width is not None:
            raise ValueError("Only provide two of the three parameters - `start`, `end` and `width` but not all!")

        if width is not None:
            if start is None and end is None:
                raise ValueError("If width is provided, either start or end must be provided.")

        narrow_ir = self._ranges.narrow(start=start, end=end, width=width)
        output = self._define_output(in_place)
        output._ranges = narrow_ir
        return output

    #####################################
    ######>> inter-range methods <<######
    #####################################

    def _group_indices_by_chrm(self, ignore_strand: bool = False) -> dict:
        grouped_indices = defaultdict(list)
        seqnames_list = self.get_seqnames(as_type="list")

        if ignore_strand:
            strands_list = ["*"] * len(self)
        else:
            strands_list = self.get_strand(as_type="list")

        for i, (seq, strand) in enumerate(zip(seqnames_list, strands_list)):
            grouped_indices[(seq, strand)].append(i)

        return dict(grouped_indices)

    def reduce(
        self,
        with_reverse_map: bool = False,
        drop_empty_ranges: bool = False,
        min_gap_width: int = 1,
        ignore_strand: bool = False,
    ) -> GenomicRanges:
        """Reduce orders the ranges, then merges overlapping or adjacent ranges.

        Args:
            with_reverse_map:
                Whether to return map of indices back to
                original object. Defaults to False.

            drop_empty_ranges:
                Whether to drop empty ranges. Defaults to False.

            min_gap_width:
                Ranges separated by a gap of
                at least ``min_gap_width`` positions are not merged. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with reduced intervals.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        _new_mcols = self._mcols.set_column("reduceindices", range(len(self)))
        _new_self = self.set_mcols(_new_mcols, in_place=False)

        all_grp_ranges = []
        rev_map_list = []
        group_keys = []

        for grp_key, grp_indices in chrm_grps.items():
            _grp_subset = _new_self[grp_indices]
            _original_indices = _grp_subset.mcols.get_column("reduceindices")

            res_ir = _grp_subset.ranges.reduce(
                with_reverse_map=True,
                drop_empty_ranges=drop_empty_ranges,
                min_gap_width=min_gap_width,
            )

            if len(res_ir) > 0:
                group_keys.extend([grp_key] * len(res_ir))
                all_grp_ranges.append(res_ir)

                if with_reverse_map:
                    for j in res_ir.mcols.get_column("revmap"):
                        rev_map_list.append([_original_indices[x] for x in j])

        if not all_grp_ranges:
            return GenomicRanges.empty()

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)
        all_merged_ranges.mcols.remove_column("revmap", in_place=True)

        new_seqnames = [k[0] for k in group_keys]
        new_strand = np.asarray([STRAND_MAP[k[1]] for k in group_keys])

        output = GenomicRanges(seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges, seqinfo=self.seqinfo)

        if with_reverse_map:
            output.mcols.set_column("revmap", rev_map_list, in_place=True)

        order = output._order_for_interranges()
        return output[order]

    def range(self, with_reverse_map: bool = False, ignore_strand: bool = False) -> GenomicRanges:
        """Calculate range bounds for each distinct (seqname, strand) pair.

        Args:
            with_reverse_map:
                Whether to return map of indices back to
                original object. Defaults to False.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with the range bounds.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        all_grp_ranges = []
        rev_map_list = []
        group_keys = []

        for grp_key, grp_indices in chrm_grps.items():
            _grp_subset = self[grp_indices]
            res_ir = _grp_subset.ranges.range()

            if len(res_ir) > 0:
                group_keys.extend([grp_key] * len(res_ir))
                all_grp_ranges.append(res_ir)
                if with_reverse_map:
                    rev_map_list.extend(grp_indices * len(res_ir))

        if not all_grp_ranges:
            return GenomicRanges.empty()

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        new_seqnames = [k[0] for k in group_keys]
        new_strand = np.asarray([STRAND_MAP[k[1]] for k in group_keys])

        output = GenomicRanges(seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges, seqinfo=self.seqinfo)

        if with_reverse_map:
            output.mcols.set_column("revmap", rev_map_list, in_place=True)

        order = output._order_for_interranges()
        return output[order]

    def gaps(
        self,
        start: int = 1,
        end: Optional[Union[int, Dict[str, int]]] = None,
        ignore_strand: bool = False,
    ) -> GenomicRanges:
        """Identify complemented ranges for each distinct (seqname, strand) pair.

        Args:
            start:
                Restrict chromosome start position. Defaults to 1.

            end:
                Restrict end position for each chromosome.
                Defaults to None. If None, extracts sequence information from
                :py:attr:`~seqinfo` object if available.

                Alternatively, you may provide a dictionary with seqnames as
                keys and the values specifying the ends.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` with complement ranges.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        all_grp_ranges = []
        group_keys = []

        all_seqs = self.seqinfo.get_seqnames()
        all_strands = ["*"] if ignore_strand else ["+", "-", "*"]

        for seq_name in all_seqs:
            for strand_str in all_strands:
                grp_key = (seq_name, strand_str)

                _end_val = None
                if isinstance(end, dict):
                    _end_val = end.get(seq_name)
                elif isinstance(end, int):
                    _end_val = end
                else:  # end is None
                    seq_idx = self.seqinfo._seqnames.index(seq_name)
                    _end_val = self.seqinfo.seqlengths[seq_idx]

                gaps_ir = None
                if grp_key in chrm_grps:
                    _grp_subset = self[chrm_grps[grp_key]]
                    gaps_ir = _grp_subset.ranges.gaps(start=start, end=_end_val)
                elif _end_val is not None:
                    # If group doesn't exist, the gap is the whole region
                    gaps_ir = IRanges([start], [_end_val - start + 1])

                if gaps_ir and len(gaps_ir) > 0:
                    all_grp_ranges.append(gaps_ir)
                    group_keys.extend([grp_key] * len(gaps_ir))

        if not all_grp_ranges:
            return GenomicRanges.empty()

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        new_seqnames = [k[0] for k in group_keys]
        new_strand = np.asarray([STRAND_MAP[k[1]] for k in group_keys])

        output = GenomicRanges(seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges, seqinfo=self.seqinfo)
        order = output._order_for_interranges()
        return output[order]

    def disjoin(self, with_reverse_map: bool = False, ignore_strand: bool = False) -> GenomicRanges:
        """Calculate disjoint genomic positions for each distinct (seqname, strand) pair.

        Args:
            with_reverse_map:
                Whether to return a map of indices back to the original object.
                Defaults to False.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` containing disjoint ranges.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        all_grp_ranges = []
        rev_map_list = []
        group_keys = []

        for grp_key, grp_indices in chrm_grps.items():
            _grp_subset = self[grp_indices]
            res_ir = _grp_subset.ranges.disjoin(with_reverse_map=True)

            if len(res_ir) > 0:
                group_keys.extend([grp_key] * len(res_ir))
                all_grp_ranges.append(res_ir)

                if with_reverse_map:
                    for j in res_ir.mcols.get_column("revmap"):
                        rev_map_list.append([grp_indices[x] for x in j])

        if not all_grp_ranges:
            return GenomicRanges.empty()

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)
        all_merged_ranges.mcols.remove_column("revmap", in_place=True)

        new_seqnames = [k[0] for k in group_keys]
        new_strand = np.asarray([STRAND_MAP[k[1]] for k in group_keys])

        output = GenomicRanges(seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges, seqinfo=self.seqinfo)

        if with_reverse_map:
            output.mcols.set_column("revmap", rev_map_list, in_place=True)

        order = output._order_for_interranges()
        return output[order]

    def is_disjoint(self, ignore_strand: bool = False) -> bool:
        """Calculate disjoint genomic positions for each distinct (seqname, strand) pair.

        Args:
            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            True if all ranges are disjoint, otherwise False.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        for grp_indices in chrm_grps.values():
            _grp_subset = self[grp_indices]
            if not _grp_subset.ranges.is_disjoint():
                return False

        return True

    def disjoint_bins(self, ignore_strand: bool = False) -> np.ndarray:
        """Split ranges into a set of bins so that the ranges in each bin are disjoint.

        Args:
            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            An ndarray indicating the bin index for each range.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        binned_results = np.zeros(len(self), dtype=int)

        for grp_indices in chrm_grps.values():
            _grp_subset = self[grp_indices]
            res_ir_bins = _grp_subset.ranges.disjoint_bins()
            binned_results[grp_indices] = res_ir_bins

        return binned_results

    def coverage(
        self, shift: int = 0, width: Optional[int] = None, weight: int = 1, ignore_strand: bool = True
    ) -> Dict[str, np.ndarray]:
        """
        Calculate coverage for each chromosome. For each position, this method
        counts the number of ranges that cover it.

        Args:
            shift:
                Shift all genomic positions. Defaults to 0.

            width:
                Restrict the width of all chromosomes. Defaults to None.
            weight:
                Weight to use for each range. Defaults to 1.

            ignore_strand:
                Whether to ignore strand, effectively combining
                all strands for coverage calculation. Defaults to True.

        Returns:
            A dictionary with chromosome names as keys and the
            coverage vector as values. The dictionary is sorted
            by chromosome name.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        sorted_keys = sorted(chrm_grps.keys())

        result = {}
        for grp_key in sorted_keys:
            grp_indices = chrm_grps[grp_key]

            cov = self._ranges[grp_indices].coverage(shift=shift, width=width, weight=weight)
            seq_name = grp_key[0]
            if seq_name in result:
                if len(result[seq_name]) < len(cov):
                    result[seq_name] = np.pad(result[seq_name], (0, len(cov) - len(result[seq_name])), "constant")
                elif len(cov) < len(result[seq_name]):
                    cov = np.pad(cov, (0, len(result[seq_name]) - len(cov)), "constant")

                result[seq_name] += cov
            else:
                result[seq_name] = cov

        return result

    def _order_for_interranges(self, decreasing: bool = False) -> np.ndarray:
        """Get the order of indices for sorting.

        Order is determined by:
        1. Chromosome name (alphabetical).
        2. Strand, in the user-specified order: '+', '-', '*'.
        3. Start position.

        Args:
            decreasing:
                Whether to sort in descending order.

        Returns:
            A NumPy vector containing index positions in the sorted order.
        """
        if len(self) == 0:
            return np.array([], dtype=int)

        strand_map = {1: 0, -1: 1, 0: 2}
        sorted_strands = np.array([strand_map[s] for s in self._strand])
        sort_keys = np.column_stack([self.start, sorted_strands, self._seqnames])
        order_indices = np.lexsort(sort_keys.T)

        if decreasing:
            return order_indices[::-1]

        return order_indices

    ################################
    ######>> set operations <<######
    ################################

    def union(self, other: GenomicRanges, ignore_strand: bool = False) -> GenomicRanges:
        """Find union of genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with all ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        grs = [self, other]
        if ignore_strand is True:
            grs = [self.set_strand(strand=None), other.set_strand(strand=None)]

        all_combined = _combine_GenomicRanges(*grs)
        output = all_combined.reduce(drop_empty_ranges=True)
        return output

    def setdiff(self, other: GenomicRanges, ignore_strand: bool = False) -> GenomicRanges:
        """Find set difference of genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with the diff ranges.
        """
        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        x = self.__copy__()
        y = other.__copy__()
        if ignore_strand is True:
            x = self.set_strand(strand=None)
            y = other.set_strand(strand=None)

        x.set_seqinfo(merge_SeqInfo((x.get_seqinfo(), y.get_seqinfo())), in_place=True)
        y.set_seqinfo(merge_SeqInfo((y.get_seqinfo(), x.get_seqinfo())), in_place=True)

        seqlengths = dict(zip(x._seqinfo.get_seqnames(), x._seqinfo.get_seqlengths()))
        all_combs = _fast_combine_GenomicRanges(x, y).range(ignore_strand=True)
        all_combs_ends = dict(zip(all_combs.get_seqnames(), all_combs.get_end()))

        for seq, seqlen in seqlengths.items():
            if seqlen is None:
                _val = all_combs_ends.get(seq, None)
                seqlengths[seq] = int(_val) if _val is not None else None

        x_gaps = x.gaps(end=seqlengths)
        x_gaps_u = x_gaps.union(y)
        diff = x_gaps_u.gaps(end=seqlengths)

        return diff

    def intersect(self, other: GenomicRanges, ignore_strand: bool = False) -> GenomicRanges:
        """Find intersecting genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object with intersecting ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        x = self.__copy__()
        y = other.__copy__()
        if ignore_strand is True:
            x = self.set_strand(strand=None)
            y = other.set_strand(strand=None)

        x.set_seqinfo(merge_SeqInfo((x.get_seqinfo(), y.get_seqinfo())), in_place=True)
        y.set_seqinfo(merge_SeqInfo((y.get_seqinfo(), x.get_seqinfo())), in_place=True)

        seqlengths = dict(zip(x._seqinfo.get_seqnames(), x._seqinfo.get_seqlengths()))
        all_combs = _fast_combine_GenomicRanges(x, y).range(ignore_strand=True)
        all_combs_ends = dict(zip(all_combs.get_seqnames(), all_combs.get_end()))

        for seq, seqlen in seqlengths.items():
            if seqlen is None:
                seqlengths[seq] = int(all_combs_ends[seq])

        y_gaps = y.gaps(end=seqlengths)
        diff = x.setdiff(y_gaps)
        return diff

    def intersect_ncls(self, other: GenomicRanges, delete_index: bool = True, num_threads: int = 1) -> GenomicRanges:
        """Find intersecting genomic intervals with `other` (uses NCLS index).

        Args:
            other:
                The other ``GenomicRanges`` object.

            delete_index:
                Defaults to True, to delete the cached ncls index.
                Set to False, to reuse the index across multiple queries.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Returns:
            A new ``GenomicRanges`` object with intersecting ranges.
        """
        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        self_end = self.get_end()
        other_end = other.get_end()

        other._ranges._build_ncls_index()
        res = other._ranges.find_overlaps(self._ranges, num_threads=num_threads)

        if delete_index:
            other._ranges.delete_nclist_index()

        _other_indexes = res["self_hits"]
        _self_indexes = res["query_hits"]
        other_chrms = np.array([other._seqinfo._seqnames[other._seqnames[i]] for i in _other_indexes])
        self_chrms = np.array([self._seqinfo._seqnames[self._seqnames[i]] for i in _self_indexes])

        other_strands = other._strand[_other_indexes]
        self_strands = self._strand[_self_indexes]

        filtered_indexes = np.logical_and(other_chrms == self_chrms, other_strands == self_strands)

        self_starts = self.start[_self_indexes][filtered_indexes]
        other_starts = other.start[_other_indexes][filtered_indexes]
        new_starts = np.where(self_starts > other_starts, self_starts, other_starts)

        self_ends = self_end[_self_indexes][filtered_indexes]
        other_ends = other_end[_other_indexes][filtered_indexes]
        new_ends = np.where(self_ends < other_ends, self_ends, other_ends)

        # filtered_keys = self_grp_keys[filtered_indexes]
        # splits = [x.split(_granges_delim) for x in filtered_keys]
        new_seqnames = self_chrms[filtered_indexes].tolist()
        new_strands = self_strands[filtered_indexes]

        output = GenomicRanges(
            seqnames=new_seqnames,
            ranges=IRanges(new_starts, new_ends - new_starts),
            strand=np.array(new_strands).astype(np.int8),
        )

        return output

    ###################################
    ######>> search operations <<######
    ###################################

    def extract_groups_by_seqnames(self):
        groups = []
        for idx, _ in enumerate(self._seqinfo._seqnames):
            idx = np.where(self._seqnames == idx)[0]
            groups.append(idx)
        return groups

    def _get_query_common_groups(self, query: GenomicRanges) -> Tuple[np.ndarray, np.ndarray]:
        # smerged = merge_SeqInfo([self._seqinfo, query._seqinfo])
        common_seqlevels = set(self._seqinfo._seqnames).intersection(query._seqinfo._seqnames)
        q_group_idx = [self._seqinfo._seqnames.index(i) for i in common_seqlevels]
        s_group_idx = [query._seqinfo._seqnames.index(i) for i in common_seqlevels]

        s_grps = self.extract_groups_by_seqnames()
        self_groups = [s_grps[i] for i in q_group_idx]

        q_grps = query.extract_groups_by_seqnames()
        query_groups = [q_grps[i] for i in s_group_idx]

        return (self_groups, query_groups)

    def find_overlaps(
        self,
        query: GenomicRanges,
        query_type: Literal["any", "start", "end", "within"] = "any",
        select: Literal["all", "first", "last", "arbitrary"] = "all",
        max_gap: int = -1,
        min_overlap: int = 0,
        ignore_strand: bool = False,
        num_threads: int = 1,
    ) -> BiocFrame:
        """Find overlaps between subject (self) and a ``query`` ``GenomicRanges`` object.

        Args:
            query:
                Query `GenomicRanges`.

            query_type:
                Overlap query type, must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".

            select:
                Determine what hit to choose when
                there are multiple hits for an interval in ``subject``.

            max_gap:
                Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).

            min_overlap:
                Minimum overlap with query. Defaults to 0.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            A BiocFrame with two columns:
            - query_hits: Indices into query ranges
            - self_hits: Corresponding indices into self ranges that are upstream
            Each row represents a query-self pair where self overlaps query.
        """
        OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        if query_type not in OVERLAP_QUERY_TYPES:
            raise ValueError(f"'{query_type}' must be one of {', '.join(OVERLAP_QUERY_TYPES)}.")

        from iranges.lib_iranges import find_overlaps_groups

        if len(self) >= len(query):
            self_groups, query_groups = self._get_query_common_groups(query)

            all_s_hits, all_q_hits = find_overlaps_groups(
                self._ranges.get_start().astype(np.int32),
                self._ranges.get_end_exclusive().astype(np.int32),
                [s.astype(np.int32) for s in self_groups],
                query._ranges.get_start().astype(np.int32),
                query._ranges.get_end_exclusive().astype(np.int32),
                [q.astype(np.int32) for q in query_groups],
                query_type,
                select,
                max_gap,
                min_overlap,
                num_threads,
            )

            if ignore_strand is False:
                s_strands = self._strand[all_s_hits]
                q_strands = query._strand[all_q_hits]

                mask = s_strands == q_strands
                # to allow '*' with any strand from query
                mask[s_strands == 0] = True
                mask[q_strands == 0] = True
                all_q_hits = all_q_hits[mask]
                all_s_hits = all_s_hits[mask]

            order = np.argsort(all_q_hits, stable=True)
            return BiocFrame({"query_hits": all_q_hits[order], "self_hits": all_s_hits[order]})
        else:
            query_groups, self_groups = query._get_query_common_groups(self)

            all_q_hits, all_s_hits = find_overlaps_groups(
                query._ranges.get_start().astype(np.int32),
                query._ranges.get_end_exclusive().astype(np.int32),
                [q.astype(np.int32) for q in query_groups],
                self._ranges.get_start().astype(np.int32),
                self._ranges.get_end_exclusive().astype(np.int32),
                [s.astype(np.int32) for s in self_groups],
                query_type,
                select,
                max_gap,
                min_overlap,
                num_threads,
            )

            if ignore_strand is False:
                q_strands = query._strand[all_q_hits]
                s_strands = self._strand[all_s_hits]

                mask = s_strands == q_strands
                mask[q_strands == 0] = True
                mask[s_strands == 0] = True
                all_q_hits = all_q_hits[mask]
                all_s_hits = all_s_hits[mask]

            order = np.argsort(all_q_hits, stable=True)
            return BiocFrame({"query_hits": all_q_hits[order], "self_hits": all_s_hits[order]})

    def count_overlaps(
        self,
        query: GenomicRanges,
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 0,
        ignore_strand: bool = False,
        num_threads: int = 1,
    ) -> np.ndarray:
        """Count overlaps between subject (self) and a ``query`` ``GenomicRanges`` object.

        Args:
            query:
                Query `GenomicRanges`.

            query_type:
                Overlap query type, must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".

            max_gap:
                Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).

            min_overlap:
                Minimum overlap with query. Defaults to 0.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            NumPy vector with length matching query,
            value represents the number of overlaps in `self` for each query.
        """
        _overlaps = self.find_overlaps(
            query,
            query_type=query_type,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
            num_threads=num_threads,
        )
        result = np.zeros(len(query))
        _ucounts = np.unique_counts(_overlaps.get_column("query_hits"))
        result[_ucounts.values] = _ucounts.counts

        return result

    def subset_by_overlaps(
        self,
        query: GenomicRanges,
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 0,
        ignore_strand: bool = False,
        num_threads: int = 1,
    ) -> GenomicRanges:
        """Subset ``subject`` (self) with overlaps in ``query`` `GenomicRanges` object.

        Args:
            query:
                Query `GenomicRanges`.

            query_type:
                Overlap query type, must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".

            max_gap:
                Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).

            min_overlap:
                Minimum overlap with query. Defaults to 0.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Raises:
            TypeError:
                If ``query`` is not of type `GenomicRanges`.

        Returns:
            A ``GenomicRanges`` object containing overlapping ranges.
        """
        _overlaps = self.find_overlaps(
            query,
            query_type=query_type,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
            num_threads=num_threads,
        )
        _all_indices = np.unique(_overlaps.get_column("self_hits"))
        return self[_all_indices]

    #################################
    ######>> nearest methods <<######
    #################################

    def nearest(
        self,
        query: GenomicRanges,
        select: Literal["all", "arbitrary"] = "arbitrary",
        ignore_strand: bool = False,
        num_threads: int = 1,
        adjacent_equals_overlap: bool = True,
    ) -> Union[np.ndarray, BiocFrame]:
        """Search nearest positions both upstream and downstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

            adjacent_equals_overlap:
                Whether to consider immediately-adjacent subject intervals to be
                equally "nearest" to the query as an overlapping subject interval.

                If true, both overlapping and immediately-adjacent subject intervals
                (i.e., a gap of zero) will be reported in matches.

                Otherwise, immediately-adjacent subjects will only be reported if
                overlapping subjects are not present.

        Returns:
            If select="arbitrary":
                A numpy array of integers with length matching query, containing indices
                into self for the closest for each query range. Value may be None if there
                are no matches.
            If select="all":
                A BiocFrame with two columns:
                - query_hits: Indices into query ranges
                - self_hits: Corresponding indices into self ranges that are upstream
                Each row represents a query-self pair where subject is nearest to query.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        from iranges.lib_iranges import nearest_groups

        self_groups, query_groups = self._get_query_common_groups(query)

        all_q_hits, all_s_hits = nearest_groups(
            self._ranges.get_start().astype(np.int32),
            self._ranges.get_end_exclusive().astype(np.int32),
            [s.astype(np.int32) for s in self_groups],
            query._ranges.get_start().astype(np.int32),
            query._ranges.get_end_exclusive().astype(np.int32),
            [q.astype(np.int32) for q in query_groups],
            select,
            num_threads,
            adjacent_equals_overlap,
        )

        if ignore_strand is False:
            s_strands = self._strand[all_s_hits]
            q_strands = query._strand[all_q_hits]

            mask = s_strands == q_strands
            # to allow '*' with any strand from query
            mask[s_strands == 0] = True
            mask[q_strands == 0] = True
            all_q_hits = all_q_hits[mask]
            all_s_hits = all_s_hits[mask]

        order = np.argsort(all_q_hits, stable=True)

        final_qhits = all_q_hits[order]
        final_shits = all_s_hits[order]

        if select == "arbitrary":
            ret_result = np.full(len(query), None)
            for s, q in zip(final_shits, final_qhits):
                ret_result[q] = s
            # unique_q, indices = np.unique(final_qhits, return_index=True)
            # ret_result[unique_q] = final_shits[indices]
            return ret_result
        else:
            return BiocFrame({"query_hits": final_qhits, "self_hits": final_shits})

    def precede(
        self,
        query: GenomicRanges,
        select: Literal["all", "first"] = "first",
        ignore_strand: bool = False,
        num_threads: int = 1,
    ) -> Union[np.ndarray, BiocFrame]:
        """Search nearest positions only downstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Whether to return "all" hits or just "first".
                Defaults to "first".

            ignore_strand:
                Whether to ignore strand. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Returns:
            If select="first":
                A numpy array of integers with length matching query, containing indices
                into self for the closest upstream position of each query range. Value may be
                None if there are no matches.
            If select="all":
                A BiocFrame with two columns:
                - query_hits: Indices into query ranges
                - self_hits: Corresponding indices into self ranges that are upstream
                Each row represents a query-self pair where self precedes query.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        effective_threads = min(num_threads, cpu_count()) if num_threads > 0 else cpu_count()
        self_groups, query_groups = self._get_query_common_groups(query)

        tasks = [
            (
                "precede",
                s_group,
                q_group,
                self._ranges,
                query._ranges,
                self._strand,
                query._strand,
                ignore_strand,
            )
            for s_group, q_group in zip(self_groups, query_groups)
        ]

        if effective_threads > 1 and len(tasks) > 1:
            with Pool(processes=effective_threads) as pool:
                results = pool.map(wrapper_follow_precede, tasks)
        else:
            results = [wrapper_follow_precede(task) for task in tasks]

        if results:
            all_qhits_list, all_shits_list = zip(*results)
        else:
            all_qhits_list, all_shits_list = [], []

        final_qhits = np.concatenate(all_qhits_list) if all_qhits_list else np.array([], dtype=np.int32)
        final_shits = np.concatenate(all_shits_list) if all_shits_list else np.array([], dtype=np.int32)

        if select == "first":
            ret_result = np.full(len(query), None)
            for s, q in zip(final_shits, final_qhits):
                if ret_result[q] is None:
                    ret_result[q] = s
            # unique_q, indices = np.unique(final_qhits, return_index=True)
            # ret_result[unique_q] = final_shits[indices]
            return ret_result
        else:
            return BiocFrame({"query_hits": final_qhits, "self_hits": final_shits})

    def follow(
        self,
        query: GenomicRanges,
        select: Literal["all", "last"] = "last",
        ignore_strand: bool = False,
        num_threads: int = 1,
    ) -> Union[np.ndarray, BiocFrame]:
        """Search nearest positions only upstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Whether to return "all" hits or just "last".
                Defaults to "last".

            ignore_strand:
                Whether to ignore strand. Defaults to False.

            num_threads:
                Number of threads to use.
                Defaults to 1.

        Returns:
            If select="last":
                A numpy array of integers with length matching query, containing indices
                into self for the closest downstream position of each query range. Value may be
                None if there are no matches.
            If select="all":
                A BiocFrame with two columns:
                - query_hits: Indices into query ranges
                - self_hits: Corresponding indices into self ranges that are upstream
                Each row represents a query-self pair where self follows query.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        effective_threads = min(num_threads, cpu_count()) if num_threads > 0 else cpu_count()
        self_groups, query_groups = self._get_query_common_groups(query)

        tasks = [
            (
                "follow",
                s_group,
                q_group,
                self._ranges,
                query._ranges,
                self._strand,
                query._strand,
                ignore_strand,
            )
            for s_group, q_group in zip(self_groups, query_groups)
        ]

        if effective_threads > 1 and len(tasks) > 1:
            with Pool(processes=effective_threads) as pool:
                results = pool.map(wrapper_follow_precede, tasks)
        else:
            results = [wrapper_follow_precede(task) for task in tasks]

        if results:
            all_qhits_list, all_shits_list = zip(*results)
        else:
            all_qhits_list, all_shits_list = [], []

        final_qhits = np.concatenate(all_qhits_list) if all_qhits_list else np.array([], dtype=np.int32)
        final_shits = np.concatenate(all_shits_list) if all_shits_list else np.array([], dtype=np.int32)

        if select == "last":
            ret_result = np.full(len(query), None)
            for q, s in zip(final_qhits, final_shits):
                ret_result[q] = s
            return ret_result
        else:
            return BiocFrame({"query_hits": final_qhits, "self_hits": final_shits})

    def distance(self, query: Union[GenomicRanges, IRanges]) -> np.ndarray:
        """Compute the pair-wise distance with intervals in query.

        Args:
            query:
                Query `GenomicRanges` or `IRanges`.

        Returns:
            Numpy vector containing distances for each interval in query.
        """
        if not isinstance(query, (IRanges, GenomicRanges)):
            raise TypeError("'query' is neither a `GenomicRanges` nor `IRanges` object.")

        if len(self) != len(query):
            raise ValueError("'query' does not contain the same number of ranges.")

        _qranges = query
        if isinstance(query, GenomicRanges):
            _qranges = query.get_ranges()

        return self._ranges.distance(_qranges)

    #############################
    ######>> comparisons <<######
    #############################

    def match(self, query: GenomicRanges, ignore_strand: bool = False) -> np.ndarray:
        """Element wise comparison to find exact match ranges.

        Args:
            query:
                Query ``GenomicRanges`` to search for matches.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Returns:
            A NumPy array with length matching query
            containing the matched indices.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        result = np.full(len(query), None)
        all_overlaps = self.find_overlaps(query, ignore_strand=ignore_strand)

        ol_query = all_overlaps.get_column("query_hits")
        ol_subject = all_overlaps.get_column("self_hits")

        unique_queries = np.unique(ol_query)

        for q in unique_queries:
            q_ctx = query[q]
            matches = np.where(ol_query == q)[0]
            self_indices = ol_subject[matches]
            self_matches = self[ol_subject[matches]]

            tgt = self_indices[
                (np.array(self_matches.get_start()) == q_ctx.get_start()[0])
                & (np.array(self_matches.get_end()) == q_ctx.get_end()[0])
            ]

            if result[q] is None:
                result[q] = tgt[0]

        return result

    def _get_ranges_as_list(self) -> List[Tuple[int, int, int]]:
        """Internal method to get ranges as a list of tuples.

        Returns:
            List of tuples containing the start, end and the index.
        """
        ranges = []
        strands = self._strand.copy()
        strands[strands == 1] = 6
        strands[strands == -1] = 7
        strands[strands == 0] = 8

        for i in range(len(self)):
            ranges.append(
                (
                    self._seqnames[i],
                    strands[i],
                    self._ranges._start[i],
                    self._ranges.end[i],
                    i,
                )
            )

        return ranges

    def order(self, decreasing: bool = False) -> np.ndarray:
        """Get the order of indices for sorting.

        Order orders the genomic ranges by chromosome and strand.
        Strand is ordered by reverse first (-1), any strand (0) and
        forward strand (-1). Then by the start positions and width if
        two regions have the same start.

        Args:
            decreasing:
                Whether to sort in descending order. Defaults to False.

        Returns:
            NumPy vector containing index positions in the sorted order.
        """
        intvals = sorted(self._get_ranges_as_list(), reverse=decreasing)
        order = [o[4] for o in intvals]
        return np.asarray(order)

    def sort(self, decreasing: bool = False, in_place: bool = False) -> GenomicRanges:
        """Get the order of indices for sorting.

        Args:
            decreasing:
                Whether to sort in descending order. Defaults to False.

            in_place:
                Whether to modify the object in place. Defaults to False.

        Returns:
            A modified ``GenomicRanges`` object with the trimmed regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        order = self.order(decreasing=decreasing)
        output = self._define_output(in_place)
        return output[list(order)]

    def rank(self) -> List[int]:
        """Get rank of the ``GenomicRanges`` object.

        For each range identifies its position is a sorted order.

        Returns:
            Numpy vector containing rank.
        """
        intvals = sorted(self._get_ranges_as_list())
        order = [o[4] for o in intvals]
        rank = [order.index(x) for x in range(len(order))]
        return rank

    ##############################
    ######>> misc methods <<######
    ##############################

    def sample(self, k: int = 5) -> GenomicRanges:
        """Randomly sample ``k`` intervals.

        Args:
            k:
                Number of intervals to sample. Defaults to 5.

        Returns:
            A new ``GenomicRanges`` with randomly sampled ranges.
        """
        from random import sample

        sample = sample(range(len(self)), k=k)
        return self[sample]

    def invert_strand(self, in_place: bool = False) -> GenomicRanges:
        """Invert strand for each range.

        Conversion map:
            - "+" map to "-"
            - "-" becomes "+"
            - "*" stays the same

        Args:
            in_place:
                Whether to modify the object in place. Defaults to False.

        Returns:
            A modified ``GenomicRanges`` object with the trimmed regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        convertor = {"1": "-1", "-1": "1", "0": "0"}
        inverts = [convertor[str(idx)] for idx in self._strand]

        output = self._define_output(in_place)
        output._strand = inverts
        return output

    ################################
    ######>> window methods <<######
    ################################

    def tile(self, n: Optional[int] = None, width: Optional[int] = None) -> List[GenomicRanges]:
        """Split each interval by ``n`` (number of sub intervals) or ``width`` (intervals with equal width).

        Note: Either ``n`` or ``width`` must be provided but not both.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a genome into chunks.

        Args:
            n:
                Number of intervals to split into.
                Defaults to None.

            width:
                Width of each interval. Defaults to None.

        Raises:
            ValueError:
                If both ``n`` and ``width`` are provided.

        Returns:
            List of ``GenomicRanges`` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either `n` or `width` must be provided but not both")

        range_tiles = self._ranges.tile(n=n, width=width)
        seqnames = self.get_seqnames("list")

        result = []
        for i in range(len(self)):
            num_tiles = len(range_tiles[i])

            result.append(
                GenomicRanges(
                    seqnames=[seqnames[i]] * num_tiles,
                    ranges=range_tiles[i],
                    strand=np.repeat(self._strand[i], num_tiles),
                    names=self._names[i] * num_tiles if self._names is not None else None,
                    seqinfo=self._seqinfo,
                )
            )

        return result

    def sliding_windows(self, width: int, step: int = 1) -> List[GenomicRanges]:
        """Slide along each range by ``width`` (intervals with equal ``width``) and ``step``.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a gneomic into chunks, or
        :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.tile_by_range`.

        Args:
            n:
                Number of intervals to split into.
                Defaults to None.

            width:
                Width of each interval. Defaults to None.

        Returns:
            A new ``GenomicRanges`` with the sliding ranges.
        """
        range_windows = self._ranges.sliding_windows(width=width, step=step)
        seqnames = self.get_seqnames("list")

        result = []
        for i in range(len(self)):
            num_tiles = len(range_windows[i])

            result.append(
                GenomicRanges(
                    seqnames=[seqnames[i]] * num_tiles,
                    ranges=range_windows[i],
                    strand=np.repeat(self._strand[i], num_tiles),
                    names=self._names[i] * num_tiles if self._names is not None else None,
                    seqinfo=self._seqinfo,
                )
            )

        return result

    # TODO: should really be a genomicrangeslist
    @classmethod
    def tile_genome(
        cls,
        seqlengths: Dict[str, int],
        ntile: Optional[int] = None,
        tilewidth: Optional[int] = None,
        cut_last_tile_in_chrom: bool = False,
    ) -> GenomicRanges:
        """Tile genome into approximately equal-sized regions.

        Args:
            seqlengths:
                Dictionary of sequence lengths by chromosome.

            ntile:
                Number of tiles (exclusive with tilewidth).

            tilewidth:
                Width of tiles (exclusive with ntile).

            cut_last_tile_in_chrom:
                Whether to cut the last tile in each chromosome.

        Returns:
            `GenomicRanges` object with ranges and bin numbers in mcols.
        """

        if ntile is not None and tilewidth is not None:
            raise ValueError("only one of 'ntile' and 'tilewidth' can be specified")

        seqlen_ = seqlengths
        if isinstance(seqlengths, SeqInfo):
            seqlen_ = dict(zip(seqlengths.get_seqnames(), seqlengths.seqlengths))

        if not seqlen_:
            raise ValueError("seqlengths must be non-empty")

        if any(length <= 0 for length in seqlen_.values()):
            raise ValueError("seqlengths contains zero or negative values")

        chroms = list(seqlen_.keys())
        lengths = np.array([seqlen_[c] for c in chroms], dtype=np.int32)
        chrom_ends = np.cumsum(lengths)
        # chrom_starts = np.r_[0, chrom_ends[:-1]]
        genome_size = chrom_ends[-1]

        if ntile is not None:
            if cut_last_tile_in_chrom:
                raise ValueError("cut_last_tile_in_chrom must be FALSE when ntile is supplied")

            if not isinstance(ntile, (int, np.integer)) or ntile < 1 or ntile > genome_size:
                raise ValueError("ntile must be an integer >= 1 and <= genome size")

            tilewidth = genome_size / ntile
            genome_breaks = np.floor(tilewidth * np.arange(1, ntile + 1)).astype(np.int32)

            all_chroms = []
            all_starts = []
            all_widths = []
            all_abs_starts = []  # track absolute starts
            all_abs_ends = []  # track absolute ends

            prev_end = 0
            for i, length in enumerate(lengths):
                curr_end = prev_end + length

                chr_mask = (genome_breaks > prev_end) & (genome_breaks <= curr_end)
                chr_breaks = genome_breaks[chr_mask]

                if len(chr_breaks) > 0:
                    # Handle first tile in chromosome
                    if not all_starts or genome_breaks[chr_mask][0] > prev_end + 1:
                        all_chroms.append(chroms[i])
                        start = 1
                        end = chr_breaks[0] - prev_end
                        all_starts.append(start)
                        all_widths.append(end)
                        all_abs_starts.append(prev_end + 1)
                        all_abs_ends.append(chr_breaks[0])

                    # Handle remaining tiles in chromosome
                    for j in range(1, len(chr_breaks)):
                        all_chroms.append(chroms[i])
                        start = chr_breaks[j - 1] - prev_end + 1
                        end = chr_breaks[j] - prev_end
                        all_starts.append(start)
                        all_widths.append(end - start + 1)
                        all_abs_starts.append(chr_breaks[j - 1] + 1)
                        all_abs_ends.append(chr_breaks[j])

                prev_end = curr_end

                if i < len(lengths) - 1 and (len(chr_breaks) == 0 or chr_breaks[-1] < curr_end):
                    all_chroms.append(chroms[i])
                    if len(chr_breaks) > 0:
                        start = chr_breaks[-1] - prev_end + length + 1
                    else:
                        start = 1
                    all_starts.append(start)
                    width = length - start + 1
                    all_widths.append(width)
                    if len(chr_breaks) > 0:
                        all_abs_starts.append(chr_breaks[-1] + 1)
                    else:
                        all_abs_starts.append(prev_end + 1)
                    all_abs_ends.append(curr_end)

            all_abs_starts = np.array(all_abs_starts)
            all_abs_ends = np.array(all_abs_ends)

            bin_width = genome_size / ntile
            bin_starts = np.arange(0, genome_size, bin_width)
            bin_ends = np.append(bin_starts[1:], genome_size)

            all_bins = np.zeros(len(all_abs_starts), dtype=np.int32)

            for i in range(len(all_abs_starts)):
                range_start = all_abs_starts[i]
                range_end = all_abs_ends[i]
                overlaps = np.minimum(bin_ends, range_end) - np.maximum(bin_starts, range_start)
                overlaps = np.maximum(overlaps, 0)
                all_bins[i] = np.argmax(overlaps) + 1
        else:
            if not isinstance(tilewidth, (int, np.integer)) or tilewidth < 1 or tilewidth > genome_size:
                raise ValueError("tilewidth must be an integer >= 1 and <= genome size")

            ntile = int(np.ceil(genome_size / tilewidth))
            return GenomicRanges.tile_genome(seqlen_, ntile=ntile, cut_last_tile_in_chrom=False)

        ranges = IRanges(start=np.array(all_starts, dtype=np.int32), width=np.array(all_widths, dtype=np.int32))
        mcols = BiocFrame({"bin": np.array(all_bins, dtype=np.int32)})
        return GenomicRanges(seqnames=all_chroms, ranges=ranges, strand=["*"] * len(all_starts), mcols=mcols)

    def binned_average(
        self,
        scorename: str,
        bins: GenomicRanges,
        outname: str = "binned_average",
        in_place: bool = False,
    ) -> GenomicRanges:
        """Calculate average for a column across all regions in ``bins``, then set a column specified by 'outname' with
        those values.

        Args:
            scorename:
                Score column to compute averages on.

            bins:
                Bins you want to use.

            outname:
                New column name to add to the object.

            in_place:
                Whether to modify ``bins`` in place.

        Raises:
            ValueError:
                If ``scorename`` column does not exist.
                ``scorename`` is not all ints or floats.
            TypeError:
                If ``bins`` is not of type `GenomicRanges`.

        Returns:
            A modified ``bins`` object with the computed averages,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """
        import statistics

        if not isinstance(bins, GenomicRanges):
            raise TypeError("'bins' is not a `GenomicRanges` object.")

        if scorename not in self._mcols.column_names:
            raise ValueError(f"'{scorename}' is not a valid column name")

        values = self._mcols.get_column(scorename)

        if not all(isinstance(x, (int, float)) for x in values):
            raise ValueError(f"'{scorename}' values must be either `ints` or `floats`.")

        outvec = []
        for _, val in bins:
            overlap = self.subset_by_overlaps(query=val)
            outvec.append(statistics.mean(overlap._mcols.get_column(scorename)))

        output = bins._define_output(in_place=in_place)
        output._mcols.set_column(outname, outvec, in_place=True)
        return output

    #######################
    ######>> split <<######
    #######################

    def split(self, groups: list) -> "CompressedGenomicRangesList":
        """Split the `GenomicRanges` object into a :py:class:`~genomicranges.grangeslist.CompressedGenomicRangesList`.

        Args:
            groups:
                A list specifying the groups or factors to split by.

                Must have the same length as the number of genomic elements
                in the object.

        Returns:
            A `CompressedGenomicRangesList` containing the groups and their
            corresponding elements.
        """

        if len(groups) != len(self):
            raise ValueError("Number of groups must match the number of genomic elements.")

        gdict = group_by_indices(groups=groups)

        _names = []
        _grs = []

        for k, v in gdict.items():
            _names.append(k)
            _grs.append(self[v])

        from .grangeslist import CompressedGenomicRangesList

        return CompressedGenomicRangesList.from_list(lst=_grs, names=_names)

    #######################
    ######>> empty <<######
    #######################

    @classmethod
    def empty(cls):
        """Create an zero-length `GenomicRanges` object.

        Returns:
            same type as caller, in this case a `GenomicRanges`.
        """
        return cls([], IRanges.empty())

    ##########################
    ######>> subtract <<######
    ##########################

    def subtract(
        self, other: GenomicRanges, min_overlap: int = 1, ignore_strand: bool = False
    ) -> "CompressedGenomicRangesList":
        """Subtract searches for features in ``x`` that overlap ``self`` by at least the number of base pairs given by
        ``min_overlap``.

        Args:
            other:
                Object to subtract.

            min_overlap:
                Minimum overlap with query.
                Defaults to 1.

            ignore_strand:
                Whether to ignore strands.
                Defaults to False.

        Returns:
            A `CompressedGenomicRangesList` with the same size as ``self`` containing
            the subtracted regions.
        """

        _other_reduce = other.reduce(ignore_strand=ignore_strand)
        hits = _other_reduce.find_overlaps(self, min_overlap=min_overlap, ignore_strand=ignore_strand)

        mapper = {}
        for i in range(len(self)):
            mapper[i] = []

        for i in range(len(hits)):
            s_hit = int(hits.get_column("self_hits")[i])
            q_hit = int(hits.get_column("query_hits")[i])

            mapper[q_hit].append(s_hit)

        mapper_with_y = {}
        for idx, val in mapper.items():
            mapper_with_y[idx] = other[val]

        psetdiff = {}
        for idx, val in mapper_with_y.items():
            if len(val) == 0:
                psetdiff[idx] = self[idx]
            else:
                psetdiff[idx] = self[idx].setdiff(val)

        from .grangeslist import CompressedGenomicRangesList

        return CompressedGenomicRangesList.from_list(lst=psetdiff.values(), names=list(psetdiff.keys()))

    ##########################
    ######>> pairwise <<######
    ##########################

    def pintersect(self, other: GenomicRanges, ignore_strand: bool = False) -> GenomicRanges:
        """Parallel intersection of genomic ranges.

        Computes the intersection for each parallel pair of ranges in ``self`` and ``other``.
        If seqnames mismatch or strands are incompatible (and not ignored), the result
        for that index is an empty range (width 0).

        Args:
            other:
                The other ``GenomicRanges`` object. Must have the same length as ``self``.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` object.
        """
        if len(self) != len(other):
            raise ValueError("'self' and 'other' must have the same length.")

        merged_seqinfo = merge_SeqInfo([self.seqinfo, other.seqinfo])

        s_names = self.get_seqnames(as_type="list")
        o_names = other.get_seqnames(as_type="list")

        s_strand = self.get_strand(as_type="numpy")
        o_strand = other.get_strand(as_type="numpy")

        new_starts = np.maximum(self.start, other.start)
        new_ends = np.minimum(self.end, other.end)

        match_seqnames = np.array([x == y for x, y in zip(s_names, o_names)])

        if not ignore_strand:
            match_strands = (s_strand * o_strand) != -1
            mask = match_seqnames & match_strands
        else:
            mask = match_seqnames

        no_overlap = new_starts > new_ends
        invalid = (~mask) | no_overlap

        final_starts = new_starts.copy()
        final_ends = new_ends.copy()

        final_starts[invalid] = 1
        final_ends[invalid] = 0

        final_widths = final_ends - final_starts + 1
        final_widths[final_widths < 0] = 0

        if ignore_strand:
            new_strands = np.zeros(len(self), dtype=int)
        else:
            new_strands = s_strand.copy()
            use_other = s_strand == 0
            new_strands[use_other] = o_strand[use_other]
            new_strands[invalid] = 0

        new_ranges = IRanges(final_starts, final_widths)

        return GenomicRanges(
            seqnames=s_names,
            ranges=new_ranges,
            strand=new_strands,
            seqinfo=merged_seqinfo,
        )


def _fast_combine_GenomicRanges(*x: GenomicRanges) -> GenomicRanges:
    return GenomicRanges(
        ranges=ut.combine_sequences(*[y._ranges for y in x]),
        seqnames=ut.combine_sequences(*[y.get_seqnames() for y in x]),
        strand=ut.combine_sequences(*[y._strand for y in x]),
        names=None,
        mcols=None,
        seqinfo=merge_SeqInfo([y._seqinfo for y in x]),
        metadata=None,
        _validate=False,
    )


@ut.combine_sequences.register(GenomicRanges)
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
        seqnames=ut.combine_sequences(*[y.get_seqnames() for y in x]),
        strand=ut.combine_sequences(*[y._strand for y in x]),
        names=all_names,
        mcols=ut.relaxed_combine_rows(*[y._mcols for y in x]),
        seqinfo=merge_SeqInfo([y._seqinfo for y in x]),
        metadata=x[0]._metadata,
        _validate=False,
    )


@ut.extract_row_names.register(GenomicRanges)
def _rownames_gr(x: GenomicRanges):
    return x.get_names()


class GRanges(GenomicRanges):
    pass
