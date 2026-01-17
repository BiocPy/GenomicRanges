from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Union

import biocutils as ut
import numpy as np
from compressed_lists import CompressedCharacterList, CompressedList, CompressedNumpyList, Partitioning
from compressed_lists.split_generic import _generic_register_helper, splitAsCompressedList
from iranges import CompressedIRangesList, IRanges

from .GenomicRanges import GenomicRanges, _combine_GenomicRanges
from .sequence_info import SeqInfo

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class CompressedGenomicRangesList(CompressedList):
    """Just as it sounds, a `CompressedGenomicRangesList` is a named-list like object.

    If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple
    genomic elements, usually where the genes start and end. Genes are themselves made of many sub
    regions, e.g. exons. `CompressedGenomicRangesList` allows us to represent this nested structure.

    Currently, this class is limited in functionality, purely a read-only class with basic accessors.

    Typical usage:

    To construct a **CompressedGenomicRangesList** object, simply pass in a list of
    :py:class:`genomicranges.GenomicRanges.GenomicRanges` objects and Optionally ``names``.

    .. code-block:: python

        a = GenomicRanges(
            seqnames=[
                "chr1",
                "chr2",
                "chr1",
                "chr3",
            ],
            ranges=IRanges(
                [
                    1,
                    3,
                    2,
                    4,
                ],
                [
                    10,
                    30,
                    50,
                    60,
                ],
            ),
            strand=[
                "-",
                "+",
                "*",
                "+",
            ],
            mcols=BiocFrame(
                {
                    "score": [
                        1,
                        2,
                        3,
                        4,
                    ]
                }
            ),
        )

        b = GenomicRanges(
            seqnames=[
                "chr2",
                "chr4",
                "chr5",
            ],
            ranges=IRanges(
                [3, 6, 4],
                [
                    30,
                    50,
                    60,
                ],
            ),
            strand=[
                "-",
                "+",
                "*",
            ],
            mcols=BiocFrame(
                {
                    "score": [
                        2,
                        3,
                        4,
                    ]
                }
            ),
        )

        grl = CompressedGenomicRangesList.from_list(
            lst=[
                gr1,
                gr2,
            ],
            names=[
                "gene1",
                "gene2",
            ],
        )

    Additionally, you may also provide metadata about the genomic elements in the dictionary
    using mcols attribute.
    """

    def __init__(
        self,
        unlist_data: GenomicRanges,
        partitioning: Partitioning,
        element_metadata: Optional[dict] = None,
        metadata: Optional[Union[Dict[str, Any], ut.NamedList]] = None,
        **kwargs,
    ):
        """Initialize a CompressedIRangesList.

        Args:
            unlist_data:
                GenomicRanges object.

            partitioning:
                Partitioning object defining element boundaries.

            element_metadata:
                Optional metadata for elements.

            metadata:
                Optional general metadata.

            kwargs:
                Additional arguments.
        """
        if not isinstance(unlist_data, GenomicRanges):
            raise TypeError("'unlist_data' is not a `GenomicRanges` object.")

        super().__init__(
            unlist_data, partitioning, element_type=GenomicRanges, element_metadata=element_metadata, metadata=metadata
        )

    @classmethod
    def from_list(
        cls,
        lst: List[GenomicRanges],
        names: Optional[Union[ut.Names, Sequence[str]]] = None,
        metadata: Optional[Union[Dict[str, Any], ut.NamedList]] = None,
    ) -> CompressedGenomicRangesList:
        """Create a `CompressedIRangesList` from a regular list.

        This concatenates the list of `GenomicRanges` objects.

        Args:
            lst:
                List of `GenomicRanges` objects.

                Must have the same number and names of columns.

            names:
                Optional names for list elements.

            metadata:
                Optional metadata.

        Returns:
            A new `CompressedList`.
        """
        unlist_data = _combine_GenomicRanges(*lst)
        partitioning = Partitioning.from_list(lst, names)
        return cls(unlist_data, partitioning, metadata=metadata)

    def extract_range(self, start: int, end: int) -> GenomicRanges:
        """Extract a range from `unlist_data`.

        This method must be implemented by subclasses to handle
        type-specific extraction from `unlist_data`.

        Args:
            start:
                Start index (inclusive).

            end:
                End index (exclusive).

        Returns:
            Extracted element.
        """
        try:
            return ut.subset_sequence(self._unlist_data, range(start, end))
        except Exception as e:
            raise NotImplementedError(
                "Custom classes should implement their own `extract_range` method for slice operations"
            ) from e

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation.
        """
        output = f"{type(self).__name__}(number_of_elements={len(self)}"
        output += ", unlist_data=" + self._unlist_data.__repr__()
        output += ", partitioning=" + self._partitioning.__repr__()
        output += (
            ", element_type=" + self._element_type.__name__
            if not isinstance(self._element_type, str)
            else self._element_type
        )

        output += ", element_metadata=" + self._element_metadata.__repr__()

        if len(self._metadata) > 0:
            output += ", metadata=" + ut.print_truncated_dict(self._metadata)

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this object.
        """
        output = f"class: {type(self).__name__}\n"

        output += f"number of elements: ({len(self)}) of type: {self._element_type.__name__ if not isinstance(self._element_type, str) else self._element_type}\n"

        output += f"unlist_data: {self._unlist_data.__str__()}\n"

        output += f"partitioning: {ut.print_truncated_list(self._partitioning)}\n"

        output += f"element_metadata({str(len(self._element_metadata))} rows): {ut.print_truncated_list(list(self._element_metadata.get_column_names()), sep=' ', include_brackets=False, transform=lambda y: y)}\n"
        output += f"metadata({str(len(self._metadata))}): {ut.print_truncated_list(list(self._metadata.keys()), sep=' ', include_brackets=False, transform=lambda y: y)}\n"

        return output

    @property
    def seqnames(self) -> CompressedCharacterList:
        """Access sequence names for each element in the list.

        Returns:
            A CompressedList containing sequence names.
        """
        return CompressedCharacterList(self._unlist_data.get_seqnames(), self._partitioning, element_type=str)

    @property
    def strand(self) -> CompressedCharacterList:
        """Access strand information for each element in the list.

        Returns:
            A CompressedList containing strand information.
        """
        return CompressedCharacterList(
            self._unlist_data.get_strand(as_type="list"),
            self._partitioning,
            element_type=str,
        )

    @property
    def ranges(self) -> CompressedIRangesList:
        """Access range information for each element in the list.

        Returns:
            A CompressedList containing IRanges objects.
        """
        return CompressedIRangesList(self._unlist_data.get_ranges(), self._partitioning, element_type=IRanges)

    @property
    def start(self) -> CompressedNumpyList:
        """Access start positions for all elements.

        Returns:
            A CompressedList of start positions.
        """
        return CompressedNumpyList(self._unlist_data.start, self._partitioning, element_type=int)

    @property
    def end(self) -> CompressedNumpyList:
        """Access end positions for all elements.

        Returns:
            A CompressedList of end positions.
        """
        return CompressedNumpyList(self._unlist_data.end, self._partitioning, element_type=int)

    @property
    def width(self) -> CompressedNumpyList:
        """Access widths for all elements.

        Returns:
            A CompressedList of widths.
        """
        return CompressedNumpyList(self._unlist_data.width, self._partitioning, element_type=int)

    @property
    def seqinfo(self) -> SeqInfo:
        """Access sequence information.

        Returns:
            The SeqInfo object from the underlying GenomicRanges.
        """
        return self._unlist_data.seqinfo

    def stack(self, index_column_name: str = "index") -> GenomicRanges:
        """Stack all elements into a single GenomicRanges object.

        Args:
            index_column_name:
                Name of the column to store the original index/names.

        Returns:
            A GenomicRanges object.
        """
        output = self._unlist_data.copy()

        # Create an index column to map back to the list element
        if self.names is not None:
            indices = np.repeat(self.names, self._partitioning.get_element_lengths())
        else:
            indices = np.repeat(range(len(self)), self._partitioning.get_element_lengths())

        new_mcols = output.mcols.set_column(index_column_name, indices)
        output.set_mcols(new_mcols, in_place=True)
        return output

    @classmethod
    def empty(cls, n: int):
        """Create an zero-length `CompressedGenomicRangesList` object.

        Args:
            n:
                Number of elements.

        Returns:
            same type as caller, in this case a `CompressedGenomicRangesList`.
        """
        _range_lengths = [0] * n

        return CompressedGenomicRangesList(
            unlist_data=GenomicRanges.empty(), partitioning=Partitioning(ends=_range_lengths)
        )


@splitAsCompressedList.register
def _(
    data: GenomicRanges,
    groups_or_partitions: Union[list, Partitioning],
    names: Optional[Union[ut.Names, Sequence[str]]] = None,
    metadata: Optional[dict] = None,
) -> CompressedGenomicRangesList:
    """Handle lists of IRanges objects."""

    partitioned_data, groups_or_partitions = _generic_register_helper(
        data=data, groups_or_partitions=groups_or_partitions, names=names
    )

    if not isinstance(partitioned_data, GenomicRanges) and len(partitioned_data) != 0:
        partitioned_data = _combine_GenomicRanges(*partitioned_data)

    return CompressedGenomicRangesList(
        unlist_data=partitioned_data, partitioning=groups_or_partitions, metadata=metadata
    )


class CompressedGRangesList(CompressedGenomicRangesList):
    pass
