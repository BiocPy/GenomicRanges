from typing import Optional, Sequence, Tuple, Union

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
    if len(strand) != num_ranges:
        raise ValueError("Length of 'strand' does not match the number of ranges.")

    if not isinstance(mcols, BiocFrame):
        raise TypeError("'mcols' is not a `BiocFrame` object.")

    if mcols.shape[0] != num_ranges:
        raise ValueError("Length of 'mcols' does not match the number of ranges.")

    if len(names) != num_ranges:
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
            seqnames = np.ndarray(
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

        if validate:
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

    #################################
    ######>> Shape and stuff <<######
    #################################

    @property
    def shape(self) -> Tuple[int, int]:
        """
        Returns:
            Tuple containing the number of rows and columns in this ``BiocFrame``.
        """
        return (len(self.seqnames), self.mcols.shape[1])

    def __len__(self) -> int:
        """
        Returns:
            Number of rows.
        """
        return self.shape[0]

    def __iter__(self) -> GenomicRangesIter:
        """Iterator over rows."""
        return GenomicRangesIter(self)

    @property
    def dims(self) -> Tuple[int, int]:
        """Alias for :py:attr:`~shape`."""
        return self.shape
