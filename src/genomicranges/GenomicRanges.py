from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union
from warnings import warn

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from .SeqInfo import SeqInfo, merge_SeqInfo
from .utils import (
    create_np_vector,
    sanitize_strand_vector,
    slide_intervals,
    split_intervals,
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

    def __init__(self, obj: "GenomicRanges") -> None:
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
            seqnames = np.asarray(
                [self._reverse_seqindex[x] for x in seqnames], dtype=np.int8
            )

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
        return current_class_const(
            seqnames=self._seqnames,
            ranges=self._ranges,
            strand=self._strand,
            names=self._names,
            mcols=self._mcols,
            seqinfo=self._seqinfo,
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

            columns = []

            header = ["seqnames", "<str>"]
            _seqnames = []
            for x in self._seqnames[indices]:
                _seqnames.append(self._seqinfo.seqnames[x])

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
                    showed = ut.truncate_strings(
                        showed, width=max(40, len(header[0]), len(header[1]))
                    )
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
    ) -> Union[Tuple[np.ndarray, List[str]], List[str]]:
        """Access sequence names.

        Args:
            as_type:
                Access seqnames as factor tuple, in which case, levels and codes
                are returned.

                If ``list``, then codes are mapped to levels and returned.

        Returns:
            List of sequence names.
        """

        if as_type == "factor":
            return self._seqnames, self._seqinfo.seqnames
        elif as_type == "list":
            return [self._seqinfo.seqnames[x] for x in self._seqnames]
        else:
            raise ValueError("Argument 'as_type' must be 'factor' or 'list'.")

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
            seqnames = np.asarray(
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

                If ``factor``, a tuple width levels as a dictionary and
                  indices to ``seqinfo.seqnames`` is returned.

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
        if len(self) == 0:
            return GenomicRanges.empty()

        idx, _ = ut.normalize_subscript(subset, len(self), self._names)

        current_class_const = type(self)
        return current_class_const(
            seqnames=self._seqnames[idx],
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
        _rdf["seqnames"] = self.get_seqnames()
        _rdf["strand"] = self.get_strand(as_type="list")

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
                Input data. must contain columns 'seqnames', 'starts' and 'widths' or "ends".

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
            width = input["ends"] - input["starts"]

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

        return cls(
            ranges=ranges, seqnames=seqnames, strand=strand, names=names, mcols=mcols
        )

    ################################
    ######>> polars interop <<######
    ################################

    def to_polars(self) -> "polars.DataFrame":
        """Convert this ``GenomicRanges`` object into a :py:class:`~polars.DataFrame`.

        Returns:
            A :py:class:`~polars.DataFrame` object.
        """
        import polars as pl

        _rdf = self._ranges.to_polars()
        _rdf = _rdf.with_columns(
            seqnames=self.get_seqnames(), strand=self.get_strand(as_type="list")
        )

        if self._names is not None:
            _rdf = _rdf.with_columns(rownames=self._names)

        if self._mcols is not None:
            if self._mcols.shape[1] > 0:
                _rdf = pl.concat([_rdf, self._mcols.to_polars()], how="horizontal")

        return _rdf

    @classmethod
    def from_polars(cls, input: "polars.DataFrame") -> "GenomicRanges":
        """Create a ``GenomicRanges`` from a :py:class:`~polars.DataFrame` object.

        Args:
            input:
                Input polars DataFrame.
                must contain columns 'seqnames', 'starts' and 'widths' or "ends".

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
            width = input["ends"] - input["starts"]

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
        mcols_df = input.drop(columns=drops)

        mcols = None
        if (not mcols_df.is_empty()) or len(mcols_df.columns) > 0:
            mcols = BiocFrame.from_polars(mcols_df)

        names = None

        return cls(
            ranges=ranges, seqnames=seqnames, strand=strand, names=names, mcols=mcols
        )

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
        width: Union[int, List[int], np.ndarray],
        fix: Literal["start", "end", "center"] = "start",
        ignore_strand: bool = False,
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Resize ranges to the specified ``width`` where either the ``start``, ``end``, or ``center`` is used as an
        anchor.

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
        _REV_FIX = {"start": "end", "end": "start", "center": "center"}
        if ignore_strand is False:
            fix = [fix] * len(self)
            for i in range(len(self.strand)):
                if self._strand[i] == -1:
                    fix[i] = _REV_FIX[fix[i]]

        output = self._define_output(in_place)
        output._ranges = self._ranges.resize(width=width, fix=fix)
        return output

    def shift(
        self, shift: Union[int, List[int], np.ndarray] = 0, in_place: bool = False
    ) -> "GenomicRanges":
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

        if shift == 0:
            return self

        output._ranges = self._ranges.shift(shift=shift)
        return output

    def promoters(
        self, upstream: int = 2000, downstream: int = 200, in_place: bool = False
    ) -> "GenomicRanges":
        """Extend intervals to promoter regions.

        Generates promoter ranges relative to the transcription start site (TSS),
        where TSS is start(x). The promoter range is expanded around the TSS
        according to the upstream and downstream arguments. Upstream represents
        the number of nucleotides in the 5' direction and downstream the number
        in the 3' direction. The full range is defined as, (`start(x) - upstream`)
        to (`start(x) + downstream - 1`).

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
        all_starts = self.start
        all_ends = self.end
        all_strands = self.strand

        start_flags = [all_strands[i] != -1 for i in range(len(all_strands))]

        new_starts = np.asarray(
            [
                (
                    all_starts[idx] - upstream
                    if start_flags[idx]
                    else all_ends[idx] - downstream
                )
                for idx in range(len(start_flags))
            ]
        )
        new_ends = np.asarray(
            [
                (
                    all_starts[idx] + downstream
                    if start_flags[idx]
                    else all_ends[idx] + upstream
                )
                for idx in range(len(start_flags))
            ]
        )

        output = self._define_output(in_place)
        output._ranges = IRanges(start=new_starts, width=(new_ends - new_starts))
        return output

    def restrict(
        self,
        start: Optional[Union[int, List[int], np.ndarray]] = None,
        end: Optional[Union[int, List[int], np.ndarray]] = None,
        keep_all_ranges: bool = False,
        in_place: bool = False,
    ) -> "GenomicRanges":
        """Restrict ranges to a given start and end positions.

        Args:
            start:
                Start position. Defaults to None.

            end:
                End position. Defaults to None.

            keep_all_ranges:
                Whether to keep intervals that do not overlap with start and end.
                Defaults to False.

            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the restricted regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """

        restricted_ir = self._ranges.restrict(
            start=start, end=end, keep_all_ranges=True
        )
        output = self._define_output(in_place)
        output._ranges = restricted_ir

        if keep_all_ranges is True:
            restricted_ir._width = np.clip(restricted_ir.width, 0, None)
        else:
            _flt_idx = list(np.where(restricted_ir.width > -1)[0])
            output = output[_flt_idx]

        return output

    def trim(self, in_place: bool = False) -> "GenomicRanges":
        """Trim sequences outside of bounds for non-circular chromosomes.

        Args:
            in_place:
                Whether to modify the ``GenomicRanges`` object in place.

        Returns:
            A modified ``GenomicRanges`` object with the trimmed regions,
            either as a copy of the original or as a reference to the
            (in-place-modified) original.
        """

        if self.seqinfo is None:
            raise ValueError("Cannot trim ranges. `seqinfo` is not available.")

        seqlengths = self.seqinfo.seqlengths
        is_circular = self.seqinfo.is_circular

        if seqlengths is None:
            raise ValueError("Cannot trim ranges. `seqlengths` is not available.")

        if is_circular is None:
            warn("considering all sequences as non-circular...")

        all_chrs = self._seqnames
        all_ends = self.end

        new_ends = []
        for i in range(len(self)):
            _t_chr = all_chrs[i]
            _end = all_ends[i]

            if (
                is_circular is not None
                and is_circular[_t_chr] is False
                and _end > seqlengths[_t_chr]
            ):
                _end = seqlengths[_t_chr] + 1

            new_ends.append(_end)

        output = self._define_output(in_place)
        output._ranges.width = np.asarray(new_ends) - output._ranges.start
        return output

    def narrow(
        self,
        start: Optional[Union[int, List[int], np.ndarray]] = None,
        width: Optional[Union[int, List[int], np.ndarray]] = None,
        end: Optional[Union[int, List[int], np.ndarray]] = None,
        in_place: bool = False,
    ) -> "GenomicRanges":
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
            raise ValueError(
                "Only provide two of the three parameters - `start`, "
                "`end` and `width` but not all!"
            )

        if width is not None:
            if start is None and end is None:
                raise ValueError(
                    "If width is provided, either start or end must be provided."
                )

        narrow_ir = self._ranges.narrow(start=start, end=end, width=width)
        output = self._define_output(in_place)
        output._ranges = narrow_ir
        return output

    def _group_indices_by_chrm(self, ignore_strand: bool = False) -> dict:
        __strand = self._strand
        if ignore_strand:
            __strand = np.zeros(len(self), dtype=np.int8)

        _seqnames = [self._seqinfo._seqnames[i] for i in self._seqnames]
        grp_keys = np.char.add(
            np.char.add(_seqnames, f"{_granges_delim}"), __strand.astype(str)
        )
        unique_grps, inverse_indices = np.unique(grp_keys, return_inverse=True)

        chrm_grps = {
            str(grp): np.where(inverse_indices == i)[0].tolist()
            for i, grp in enumerate(unique_grps)
        }

        return chrm_grps

    def reduce(
        self,
        with_reverse_map: bool = False,
        drop_empty_ranges: bool = False,
        min_gap_width: int = 1,
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
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

        _new_mcols = self._ranges._mcols.set_column("reduceindices", range(len(self)))
        _new_ranges = self._ranges.set_mcols(_new_mcols)
        _new_self = self.set_ranges(_new_ranges)

        all_grp_ranges = []
        rev_map = []
        groups = []

        for seq in _new_self._seqinfo._seqnames:
            _iter_strands = [0] if ignore_strand is True else [1, -1, 0]
            for strd in _iter_strands:
                _key = f"{seq}{_granges_delim}{strd}"
                if _key in chrm_grps:
                    _grp_subset = _new_self[chrm_grps[_key]]
                    _oindices = _grp_subset._ranges._mcols.get_column("reduceindices")

                    res_ir = _grp_subset._ranges.reduce(
                        with_reverse_map=True,
                        drop_empty_ranges=drop_empty_ranges,
                        min_gap_width=min_gap_width,
                    )

                    groups.extend([_key] * len(res_ir))
                    all_grp_ranges.append(res_ir)
                    _rev_map = []
                    for j in res_ir._mcols.get_column("revmap"):
                        _rev_map.append([_oindices[x] for x in j])
                    rev_map.extend(_rev_map)

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        splits = [x.split(_granges_delim) for x in groups]
        new_seqnames = [x[0] for x in splits]
        new_strand = np.asarray([int(x[1]) for x in splits])

        output = GenomicRanges(
            seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges
        )

        if with_reverse_map is True:
            output._mcols.set_column("revmap", rev_map, in_place=True)

        # self._ranges._mcols.remove_column("reduceindices", in_place=True)

        return output

    def range(
        self, with_reverse_map: bool = False, ignore_strand: bool = False
    ) -> "GenomicRanges":
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
        rev_map = []
        groups = []

        for seq in self._seqinfo.seqnames:
            _iter_strands = [0] if ignore_strand is True else [1, -1, 0]
            for strd in _iter_strands:
                _key = f"{seq}{_granges_delim}{strd}"
                if _key in chrm_grps:
                    _grp_subset = self[chrm_grps[_key]]
                    res_ir = _grp_subset._ranges.range()

                    groups.extend([_key] * len(res_ir))
                    all_grp_ranges.append(res_ir)
                    rev_map.extend(chrm_grps[_key] * len(res_ir))

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        splits = [x.split(_granges_delim) for x in groups]
        new_seqnames = [x[0] for x in splits]
        new_strand = np.asarray([int(x[1]) for x in splits])

        output = GenomicRanges(
            seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges
        )

        if with_reverse_map is True:
            output._mcols.set_column("revmap", rev_map, in_place=True)

        return output

    def gaps(
        self,
        start: int = 1,
        end: Optional[Union[int, Dict[str, int]]] = None,
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
        """Identify complemented ranges for each distinct (seqname, strand) pair.

        Args:
            start:
                Restrict chromosome start position. Defaults to 1.

            end:
                Restrict end position for each chromosome.
                Defaults to None. If None, extracts sequence information from
                :py:attr:`~seqinfo` object if available.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Returns:
            A new ``GenomicRanges`` with complement ranges.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        all_grp_ranges = []
        groups = []

        for i, chrm in enumerate(self._seqinfo.seqnames):
            _iter_strands = [0] if ignore_strand is True else [1, -1, 0]
            for strd in _iter_strands:
                _key = f"{chrm}{_granges_delim}{strd}"

                _end = None
                if isinstance(end, dict):
                    _end = end[chrm]
                elif isinstance(end, int):
                    _end = end

                gaps = None
                if _key in chrm_grps:
                    _grp_subset = self[chrm_grps[_key]]
                    gaps = _grp_subset._ranges.gaps(
                        start=start,
                        end=_end,  # - 1 if _end is not None else _end
                    )
                else:
                    if _end is None:
                        _end = self._seqinfo.seqlengths[i]

                    if _end is not None:
                        gaps = IRanges(start=[start], width=[_end - start + 1])

                if gaps is not None:
                    all_grp_ranges.append(gaps)
                    groups.extend([_key] * len(gaps))

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        splits = [x.split(_granges_delim) for x in groups]
        new_seqnames = [x[0] for x in splits]
        new_strand = np.asarray([int(x[1]) for x in splits])

        output = GenomicRanges(
            seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges
        )

        return output

    def disjoin(
        self, with_reverse_map: bool = False, ignore_strand: bool = False
    ) -> "GenomicRanges":
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
        rev_map = []
        groups = []

        for seq in self._seqinfo.seqnames:
            _iter_strands = [0] if ignore_strand is True else [1, -1, 0]
            for strd in _iter_strands:
                _key = f"{seq}{_granges_delim}{strd}"
                if _key in chrm_grps:
                    _grp_subset = self[chrm_grps[_key]]
                    res_ir = _grp_subset._ranges.disjoin(with_reverse_map=True)

                    groups.append(_key)
                    all_grp_ranges.append(res_ir)

                    _rev_map = []
                    for j in res_ir._mcols.get_column("revmap"):
                        _rev_map.append([chrm_grps[_key][x] for x in j])
                    rev_map.append(_rev_map[0])

        all_merged_ranges = ut.combine_sequences(*all_grp_ranges)

        splits = [x.split(_granges_delim) for x in groups]
        new_seqnames = [x[0] for x in splits]
        new_strand = np.asarray([int(x[1]) for x in splits])

        output = GenomicRanges(
            seqnames=new_seqnames, strand=new_strand, ranges=all_merged_ranges
        )

        if with_reverse_map is True:
            output._mcols.set_column("revmap", rev_map, in_place=True)

        return output

    def coverage(
        self, shift: int = 0, width: Optional[int] = None, weight: int = 1
    ) -> Dict[str, np.ndarray]:
        """Calculate coverage for each chromosome, For each position, counts the number of ranges that cover it.

        Args:
            shift:
                Shift all genomic positions. Defaults to 0.

            width:
                Restrict the width of all
                chromosomes. Defaults to None.

            weight:
                Weight to use. Defaults to 1.

        Returns:
            A dictionary with chromosome names as keys and the
            coverage vector as value.
        """
        chrm_grps = self._group_indices_by_chrm(ignore_strand=True)

        shift_arr = None
        if shift > 0:
            shift_arr = np.zeros(shift)

        result = {}
        for chrm, group in chrm_grps.items():
            _grp_subset = self[group]

            all_intvals = [
                (x[0], x[1])
                for x in zip(_grp_subset._ranges._start, _grp_subset._ranges.end)
            ]

            cov, _ = create_np_vector(intervals=all_intvals, with_reverse_map=False)

            if shift > 0:
                cov = ut.combine_sequences(shift_arr, cov)

            if weight > 0:
                cov = cov * weight

            if width is not None:
                cov = cov[:width]

            result[chrm.split(_granges_delim)[0]] = cov

        return result

    ################################
    ######>> set operations <<######
    ################################

    def union(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find union of genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

        Raises:
            TypeError:
                If ``other`` is not a ``GenomicRanges``.

        Returns:
            A new ``GenomicRanges`` object with all ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        all_combined = _combine_GenomicRanges(self, other)
        output = all_combined.reduce(min_gap_width=0, drop_empty_ranges=True)
        return output

    def setdiff(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find set difference of genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

        Raises:
            TypeError:
                If ``other`` is not of type ``GenomicRanges``.


        Returns:
            A new ``GenomicRanges`` object with the diff ranges.
        """
        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        all_combined = _fast_combine_GenomicRanges(self, other)
        range_bounds = all_combined.range(ignore_strand=True)
        rb_ends = {}
        for _, val in range_bounds:
            rb_ends[val.seqnames[0]] = val.end[0]

        x_gaps = self.gaps(end=rb_ends)
        x_gaps_u = x_gaps.union(other)
        diff = x_gaps_u.gaps(end=rb_ends)

        return diff

    def intersect(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find intersecting genomic intervals with `other`.

        Args:
            other:
                The other ``GenomicRanges`` object.

        Raises:
            TypeError:
                If ``other`` is not a ``GenomicRanges``.

        Returns:
            A new ``GenomicRanges`` object with intersecting ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        all_combined = _fast_combine_GenomicRanges(self, other)
        range_bounds = all_combined.range(ignore_strand=True)
        rb_ends = {}
        for _, val in range_bounds:
            rb_ends[val.seqnames[0]] = val.end[0]

        _gaps = other.gaps(end=rb_ends)
        # _inter = self.setdiff(_gaps)
        x_gaps = self.gaps(end=rb_ends)
        x_gaps_u = x_gaps.union(_gaps)
        diff = x_gaps_u.gaps(end=rb_ends)
        return diff

    def _get_chrm_grps(self):
        chrm_grps = []
        for i in range(len(self)):
            __strand = self._strand[i]
            _grp = f"{self._seqinfo._seqnames[self._seqnames[i]]}{_granges_delim}{__strand}"
            chrm_grps.append(_grp)
        return chrm_grps

    def intersect_ncls(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find intersecting genomic intervals with `other` (uses NCLS index).

        Args:
            other:
                The other ``GenomicRanges`` object.

        Raises:
            TypeError:
                If ``other`` is not a ``GenomicRanges``.

        Returns:
            A new ``GenomicRanges`` object with intersecting ranges.
        """
        if not isinstance(other, GenomicRanges):
            raise TypeError("'other' is not a `GenomicRanges` object.")

        if not ut.package_utils.is_package_installed("ncls"):
            raise ImportError("package: 'ncls' is not installed.")

        from ncls import NCLS

        other_ncls = NCLS(other.start, other.end, np.arange(len(other)))
        _self_indexes, _other_indexes = other_ncls.all_overlaps_both(
            self.start, self.end, np.arange(len(self))
        )

        other_grp_keys = np.array(other._get_chrm_grps())
        self_grp_keys = np.array(self._get_chrm_grps())
        all_self_keys = self_grp_keys[_self_indexes]
        filtered_indexes = all_self_keys == other_grp_keys[_other_indexes]

        self_starts = self.start[_self_indexes][filtered_indexes]
        other_starts = other.start[_other_indexes][filtered_indexes]
        new_starts = np.where(self_starts > other_starts, self_starts, other_starts)

        self_ends = self.end[_self_indexes][filtered_indexes]
        other_ends = other.end[_other_indexes][filtered_indexes]
        new_ends = np.where(self_ends < other_ends, self_ends, other_ends)

        filtered_keys = all_self_keys[filtered_indexes]
        splits = [x.split(_granges_delim) for x in filtered_keys]
        new_seqnames, new_strands = zip(*splits)

        output = GenomicRanges(
            seqnames=new_seqnames,
            ranges=IRanges(new_starts, new_ends - new_starts),
            strand=np.array(new_strands).astype(np.int8),
        )

        return output

    ###################################
    ######>> search operations <<######
    ###################################

    def find_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        select: Literal["all", "first", "last", "arbitrary"] = "all",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> List[List[int]]:
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
                Minimum overlap with query. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            A list with the same length as ``query``,
            containing hits to overlapping indices.
        """
        OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        if query_type not in OVERLAP_QUERY_TYPES:
            raise ValueError(
                f"'{query_type}' must be one of {', '.join(OVERLAP_QUERY_TYPES)}."
            )

        rev_map = [[] for _ in range(len(query))]
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.find_overlaps(
                    query=_query_subset._ranges,
                    query_type=query_type,
                    select=select,
                    max_gap=max_gap,
                    min_overlap=min_overlap,
                    delete_index=False,
                )

                for j, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map[indices[j]] = _rev_map

        return rev_map

    def count_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> List[int]:
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
                Minimum overlap with query. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            A list with the same length as ``query``,
            containing number of overlapping indices.
        """
        OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        if query_type not in OVERLAP_QUERY_TYPES:
            raise ValueError(
                f"'{query_type}' must be one of {', '.join(OVERLAP_QUERY_TYPES)}."
            )

        rev_map = [0 for _ in range(len(query))]
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.find_overlaps(
                    query=_query_subset._ranges,
                    query_type=query_type,
                    select="all",
                    max_gap=max_gap,
                    min_overlap=min_overlap,
                    delete_index=False,
                )

                for j, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map[indices[j]] = len(_rev_map)

        return rev_map

    def subset_by_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
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
                Minimum overlap with query. Defaults to 1.

            ignore_strand:
                Whether to ignore strands. Defaults to False.

        Raises:
            TypeError:
                If ``query`` is not of type `GenomicRanges`.

        Returns:
            A ``GenomicRanges`` object containing overlapping ranges.
        """
        OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        if query_type not in OVERLAP_QUERY_TYPES:
            raise ValueError(
                f"'{query_type}' must be one of {', '.join(OVERLAP_QUERY_TYPES)}."
            )

        rev_map = []
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.find_overlaps(
                    query=_query_subset._ranges,
                    query_type=query_type,
                    select="all",
                    max_gap=max_gap,
                    min_overlap=min_overlap,
                    delete_index=False,
                )

                for _, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map.extend(_rev_map)

        return self[list(set(rev_map))]

    #################################
    ######>> nearest methods <<######
    #################################

    def nearest(
        self,
        query: "GenomicRanges",
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> List[List[int]]:
        """Search nearest positions both upstream and downstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Returns:
            A list with the same length as ``query``,
            containing hits to nearest indices.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        rev_map = [[] for _ in range(len(query))]
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.nearest(
                    query=_query_subset._ranges, select=select, delete_index=False
                )

                for j, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map[indices[j]] = _rev_map

        return rev_map

    def precede(
        self,
        query: "GenomicRanges",
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> List[List[int]]:
        """Search nearest positions only downstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Returns:
            A List with the same length as ``query``,
            containing hits to nearest indices.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        rev_map = [[] for _ in range(len(query))]
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.precede(
                    query=_query_subset._ranges, select=select, delete_index=False
                )

                for j, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map[indices[j]] = _rev_map

        return rev_map

    def follow(
        self,
        query: "GenomicRanges",
        select: Literal["all", "arbitrary"] = "all",
        ignore_strand: bool = False,
    ) -> List[List[int]]:
        """Search nearest positions only upstream that overlap with each range in ``query``.

        Args:
            query:
                Query ``GenomicRanges`` to find nearest positions.

            select:
                Determine what hit to choose when there are
                multiple hits for an interval in ``query``.

            ignore_strand:
                Whether to ignore strand. Defaults to False.

        Returns:
            A List with the same length as ``query``,
            containing hits to nearest indices.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        rev_map = [[] for _ in range(len(query))]
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)
        query_chrm_grps = query._group_indices_by_chrm(ignore_strand=ignore_strand)

        for group, indices in query_chrm_grps.items():
            if group in subject_chrm_grps:
                _sub_subset = self[subject_chrm_grps[group]]
                _query_subset = query[indices]

                res_idx = _sub_subset._ranges.follow(
                    query=_query_subset._ranges, select=select, delete_index=False
                )

                for j, val in enumerate(res_idx):
                    _rev_map = [subject_chrm_grps[group][x] for x in val]
                    rev_map[indices[j]] = _rev_map

        return rev_map

    def distance(self, query: Union["GenomicRanges", IRanges]) -> np.ndarray:
        """Compute the pair-wise distance with intervals in query.

        Args:
            query:
                Query `GenomicRanges` or `IRanges`.

        Returns:
            Numpy vector containing distances for each interval in query.
        """
        if not isinstance(query, (IRanges, GenomicRanges)):
            raise TypeError("'query' is not a `GenomicRanges` or `IRanges` object.")

        if len(self) != len(query):
            raise ValueError("'query' does not contain the same number of intervals.")

        _qranges = query
        if isinstance(query, GenomicRanges):
            _qranges = query.get_ranges()

        return self._ranges.distance(_qranges)

    def match(self, query: "GenomicRanges") -> List[List[int]]:
        """Element wise comparison to find exact match ranges.

        Args:
            query:
                Query ``GenomicRanges`` to search for matches.

        Raises:
            TypeError:
                If ``query`` is not of type ``GenomicRanges``.

        Returns:
            A List with the same length as ``query``,
            containing hits to matching indices.
        """
        if not isinstance(query, GenomicRanges):
            raise TypeError("'query' is not a `GenomicRanges` object.")

        ignore_strand = False
        subject_chrm_grps = self._group_indices_by_chrm(ignore_strand=ignore_strand)

        rev_map = []
        groups = []

        for i in range(len(query)):
            try:
                _seqname = query.seqnames[i]
            except Exception as _:
                warn(f"'{query.seqnames[i]}' is not present in subject.")

            _strand = query._strand[i]

            if ignore_strand is True:
                _strand = 0

            _key = f"{_seqname}{_granges_delim}{_strand}"
            if _key in subject_chrm_grps:
                _grp_subset = self[subject_chrm_grps[_key]]

                res_idx = _grp_subset._ranges.find_overlaps(
                    query=query._ranges[i],
                    query_type="any",
                    select="all",
                    delete_index=False,
                )

                groups.append(i)

                _rev_map = []
                for j in res_idx:
                    for x in j:
                        _mrange = self[subject_chrm_grps[_key][x]]._ranges

                        if (
                            _mrange.start[0] == query._ranges[i].start[0]
                            and _mrange.width[0] == query._ranges[i].width[0]
                        ):
                            _rev_map.append(subject_chrm_grps[_key][x])
                rev_map.append(_rev_map)

        return rev_map

    def _get_ranges_as_list(self) -> List[Tuple[int, int, int]]:
        """Internal method to get ranges as a list of tuples.

        Returns:
            List of tuples containing the start, end and the index.
        """
        ranges = []
        for i in range(len(self)):
            ranges.append(
                (
                    self._seqnames[i],
                    self._strand[i],
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

    def sort(self, decreasing: bool = False, in_place: bool = False) -> "GenomicRanges":
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

    def sample(self, k: int = 5) -> "GenomicRanges":
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

    def invert_strand(self, in_place: bool = False) -> "GenomicRanges":
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

    def tile_by_range(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each sequence length into chunks by ``n`` (number of intervals) or ``width`` (intervals with equal
        width).

        Note: Either ``n`` or ``width`` must be provided, but not both.

        Also, checkout :py:meth:`~genomicranges.io.tiling.tile_genome` for splitting
        the genome into chunks.

        Args:
            n:
                Number of intervals to split into.
                Defaults to None.

            width:
                Width of each interval. Defaults to None.

        Raises:
            ValueError:
                If both ``n`` or ``width`` are provided.

        Returns:
            A new ``GenomicRanges`` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("Either `n` or `width` must be provided but not both.")

        ranges = self.range()

        seqnames = []
        strand = []
        starts = []
        widths = []

        for _, val in ranges:
            twidth = None
            if n is not None:
                twidth = int((val._ranges.width[0]) / n)
            elif width is not None:
                twidth = width

            all_intervals = split_intervals(
                val._ranges._start[0], val._ranges.end[0] - 1, twidth
            )

            seqnames.extend([val.seqnames[0]] * len(all_intervals))
            strand.extend([int(val.strand[0])] * len(all_intervals))
            starts.extend([x[0] for x in all_intervals])
            widths.extend(x[1] for x in all_intervals)

        return GenomicRanges(
            seqnames=seqnames, strand=strand, ranges=IRanges(start=starts, width=widths)
        )

    def tile(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each interval by ``n`` (number of sub intervals) or ``width`` (intervals with equal width).

        Note: Either ``n`` or ``width`` must be provided but not both.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a genome into chunks, or
        :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.tile_by_range`.

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
            A new ``GenomicRanges`` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either `n` or `width` must be provided but not both")

        import math

        seqnames = []
        strand = []
        starts = []
        widths = []

        counter = 0
        for _, val in self:
            twidth = None
            if n is not None:
                twidth = math.ceil((val._ranges._width + 1) / (n))

                if twidth < 1:
                    raise RuntimeError(
                        f"'width' of region is less than 'n' for range in: {counter}."
                    )
            elif width is not None:
                twidth = width

            all_intervals = split_intervals(
                val._ranges._start[0], val._ranges.end[0] - 1, twidth
            )

            seqnames.extend([val.seqnames[0]] * len(all_intervals))
            strand.extend([int(val.strand[0])] * len(all_intervals))
            starts.extend([x[0] for x in all_intervals])
            widths.extend(x[1] for x in all_intervals)

            counter += 1

        return GenomicRanges(
            seqnames=seqnames, strand=strand, ranges=IRanges(start=starts, width=widths)
        )

    def sliding_windows(self, width: int, step: int = 1) -> "GenomicRanges":
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
        seqnames = []
        strand = []
        starts = []
        widths = []

        for _, val in self:
            all_intervals = slide_intervals(
                val._ranges._start[0],
                val._ranges.end[0] - 1,
                width=width,
                step=step,
            )

            seqnames.extend([val.seqnames[0]] * len(all_intervals))
            strand.extend([int(val.strand[0])] * len(all_intervals))
            starts.extend([x[0] for x in all_intervals])
            widths.extend(x[1] for x in all_intervals)

        return GenomicRanges(
            seqnames=seqnames, strand=strand, ranges=IRanges(start=starts, width=widths)
        )

    @classmethod
    def tile_genome(
        cls,
        seqlengths: Union[Dict, SeqInfo],
        n: Optional[int] = None,
        width: Optional[int] = None,
    ) -> "GenomicRanges":
        """Create a new ``GenomicRanges`` by partitioning a specified genome.

        If ``n`` is provided, the region is split into ``n`` intervals. The last interval may
        not contain the same 'width' as the other regions.

        Alternatively, ``width`` may be provided for each interval. Similarly, the last region
        may be less than ``width``.

        Either ``n`` or ``width`` must be provided but not both.

        Args:
            seqlengths:
                Sequence lengths of each chromosome.

                ``seqlengths`` may be a dictionary, where keys specify the chromosome
                and the value is the length of each chromosome in the genome.

                Alternatively, ``seqlengths`` may be an instance of
                :py:class:`~genomicranges.SeqInfo.SeqInfo`.

            n:
                Number of intervals to split into.
                Defaults to None, then 'width' of each interval is computed from ``seqlengths``.

            width:
                Width of each interval. Defaults to None.

        Raises:
            ValueError:
                If both ``n`` and ``width`` are provided.

        Returns:
            A new ``GenomicRanges`` with the tiled regions.
        """
        import math

        if n is not None and width is not None:
            raise ValueError("Both `n` or `width` are provided!")

        seqlen_ = seqlengths
        if isinstance(seqlengths, SeqInfo):
            seqlen_ = dict(zip(seqlengths.seqnames, seqlengths.seqlengths))

        seqnames = []
        strand = []
        starts = []
        widths = []
        for chrm, chrlen in seqlen_.items():
            twidth = None
            if n is not None:
                twidth = math.ceil(chrlen / n)
            elif width is not None:
                twidth = width

            all_intervals = split_intervals(1, chrlen, twidth)

            seqnames.extend([chrm] * len(all_intervals))
            strand.extend(["*"] * len(all_intervals))
            starts.extend([x[0] for x in all_intervals])
            widths.extend(x[1] for x in all_intervals)

        return GenomicRanges(
            seqnames=seqnames, strand=strand, ranges=IRanges(start=starts, width=widths)
        )

    def binned_average(
        self,
        scorename: str,
        bins: "GenomicRanges",
        outname: str = "binned_average",
        in_place: bool = False,
    ) -> "GenomicRanges":
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

    @classmethod
    def empty(cls):
        """Create an zero-length `GenomicRanges` object.

        Returns:
            same type as caller, in this case a `GenomicRanges`.
        """
        return cls([], IRanges.empty())


def _fast_combine_GenomicRanges(*x: GenomicRanges) -> GenomicRanges:
    return GenomicRanges(
        ranges=ut.combine_sequences(*[y._ranges for y in x]),
        seqnames=ut.combine_sequences(*[y._seqnames for y in x]),
        strand=ut.combine_sequences(*[y._strand for y in x]),
        names=None,
        mcols=None,
        seqinfo=merge_SeqInfo([y._seqinfo for y in x]),
        metadata=None,
        validate=False,
    )


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
        mcols=ut.relaxed_combine_rows(*[y._mcols for y in x]),
        seqinfo=merge_SeqInfo([y._seqinfo for y in x]),
        metadata=x[0]._metadata,
        validate=False,
    )
