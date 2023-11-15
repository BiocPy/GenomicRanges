from typing import Optional, Sequence, Union

import biocutils as ut
import numpy as np
from biocframe import BiocFrame
from iranges import IRanges

from .SeqInfo import SeqInfo
from .utils import make_strand_vector

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _guess_num_ranges(seqnames, ranges):
    if len(seqnames) != len(ranges):
        raise ValueError("Number of 'seqnames' and 'ranges' do not match!")

    return len(seqnames)


def _validate_seq_info(seq_info):
    if not isinstance(seq_info, SeqInfo):
        raise TypeError("'seq_info' is not an instance of `SeqInfo` class.")


def _validate_ranges(ranges):
    if not isinstance(ranges, IRanges):
        raise TypeError("'ranges' is not an instance of `IRanges`.")


def _validate_optional_attrs(num_ranges, strand, mcols, names):
    if strand is not None and len(strand) != num_ranges:
        raise ValueError("Length of 'strand' does not match the number of ranges.")

    if mcols is not None:
        if not isinstance(mcols, BiocFrame):
            raise TypeError("'mcols' is not a `BiocFrame` object.")

        if mcols.shape[0] != num_ranges:
            raise ValueError("Length of 'mcols' does not match the number of ranges.")

    if names is not None and len(names) != num_ranges:
        raise ValueError("Length of 'names' does not match the number of ranges.")


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

            iter_slice = self._gr.row(self._current_index)
            self._current_index += 1
            return (iter_row_index, iter_slice)

        raise StopIteration


class GenomicRanges:
    """``GenomicRanges`` provides a container class to represent and operate over genomic regions and annotations.

    Additionally, (checkout
    :py:class:`~genomicranges.SeqInfo.SeqInfo`) might also contain metadata about the
    genome, e.g. if it's circular (`is_circular`) or not.

    Note: The documentation for some of the methods are derived from the
    `GenomicRanges R/Bioconductor package <https://github.com/Bioconductor/GenomicRanges>`_.
    """

    required_columns = ["seqnames", "starts", "ends", "strand"]

    def __init__(
        self,
        seqnames: Union[Sequence[str], np.ndarray],
        ranges: IRanges,
        strand: Optional[Union[Sequence[str], Sequence[int], np.ndarray]] = None,
        names: Optional[Sequence[str]] = None,
        mcols: Optional[BiocFrame] = None,
        seq_info: Optional[SeqInfo] = None,
        metadata: Optional[dict] = None,
        validate: bool = True,
    ):
        """Initialize a ``GenomicRanges`` object.

        Args:
            seqnames:
                Sequence or chromosome names.

            ranges:
                Genomic positions and widths of the range. Must have the same length as ``seqnames``.

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

            seq_info:
                Sequence information. Defaults to None, in which case a
                :py:class:`~genomicranges.SeqInfo.SeqInfo` object is created with the unique set of
                chromosome names from ``seqnames``.

            metadata:
                Additional metadata. Defaults to None, and is assigned to an empty dictionary.

            validate:
                Internal use only.
        """
        if seq_info is None:
            seq_info = SeqInfo(seqnames=list(set(seqnames)))
        self._seq_info = seq_info

        if not isinstance(seqnames, np.ndarray):
            seqnames = np.array(
                [self._seq_info.seqnames.index(x) for x in list(seqnames)]
            )
        self._seqnames = seqnames

        self._ranges = ranges

        if strand is None:
            strand = np.repeat(0, len(self._seqnames))
        else:
            strand = make_strand_vector(strand)
        self._strand = strand

        if names is not None and not isinstance(names, ut.StringList):
            names = ut.StringList(names)
        self._names = names

        if mcols is None:
            mcols = BiocFrame({}, number_of_rows=len(self._seqnames))
        self._mcols = mcols

        self._metadata = metadata if metadata is not None else {}

        if validate is True:
            _validate_seq_info(self._seq_info)
            _validate_ranges(self._ranges)
            num_ranges = _guess_num_ranges(self._seqnames, self._ranges)
            _validate_optional_attrs(num_ranges, self._strand, self._mcols, self._names)

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
        _seqinfo_copy = deepcopy(self._seq_info)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            seqnames=_seqnames_copy,
            ranges=_ranges_copy,
            strand=_strand_copy,
            names=_names_copy,
            mcols=_mcols_copy,
            seq_info=_seqinfo_copy,
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
            seq_info=self._seq_info,
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

        if self._seq_info is not None:
            output += ", seq_info" + repr(self._seq_info)

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
            for col in self._column_names:
                data = self._data[col]
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
        if self.column_data is not None and self.column_data.shape[1]:
            footer.append(
                "column_data("
                + str(self.column_data.shape[1])
                + "): "
                + ut.print_truncated_list(
                    self.column_data.column_names,
                    sep=" ",
                    include_brackets=False,
                    transform=lambda y: y,
                )
            )
        if len(self.metadata):
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
