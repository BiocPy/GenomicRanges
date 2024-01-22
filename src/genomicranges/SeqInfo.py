from typing import Dict, List, Optional, Sequence, Union
from warnings import warn

import biocutils as ut

from .utils import _sanitize_vec

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def _validate_seqnames(seqnames):
    if not ut.is_list_of_type(seqnames, str):
        raise ValueError("'seqnames' should be a list of strings.")

    n = len(seqnames)
    if n != len(set(seqnames)):
        raise ValueError("'seqnames' should contain unique strings.")


def _validate_seqlengths(seqlengths, num_seqs):
    if not ut.is_list_of_type(seqlengths, int, ignore_none=True):
        raise ValueError("'seqlengths' should be a list of integers.")

    if num_seqs != len(seqlengths):
        raise ValueError("'seqnames' and 'seqlengths' should have the same length.")

    for sl in seqlengths:
        if sl is not None and sl < 0:
            raise ValueError("all entries of 'seqlengths' should be non-negative.")


def _validate_is_circular(is_circular, num_seqs):
    if not ut.is_list_of_type(is_circular, bool, ignore_none=True):
        raise ValueError("'is_circular' should be a list of booleans.")

    if num_seqs != len(is_circular):
        raise ValueError("'seqnames' and 'is_circular' should have the same length.")


def _validate_genome(genome, num_seqs):
    if not ut.is_list_of_type(genome, str, ignore_none=True):
        raise ValueError("'genome' should be a list of strings.")

    if num_seqs != len(genome):
        raise ValueError("'seqnames' and 'genome' should have the same length.")


class SeqInfoIterator:
    """An iterator to a :py:class:`~SeqInfo` object."""

    def __init__(self, obj: "SeqInfo") -> None:
        """Initialize the iterator.

        Args:
            obj:
                Source object to iterate.
        """
        self._sinfo = obj
        self._current_index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._current_index < len(self._sinfo):
            iter_row_index = self._sinfo.seqnames[self._current_index]

            iter_slice = self._sinfo[self._current_index]
            self._current_index += 1
            return (iter_row_index, iter_slice)

        raise StopIteration


class SeqInfo:
    """Information about the reference sequences, specifically the name and length of each sequence, whether it is a
    circular, and the identity of the genome from which it was derived."""

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Optional[Union[int, Sequence[int], Dict[str, int]]] = None,
        is_circular: Optional[Union[bool, Sequence[bool], Dict[str, bool]]] = None,
        genome: Optional[Union[str, Sequence[str], Dict[str, str]]] = None,
        validate: bool = True,
    ) -> None:
        """
        Args:
            seqnames:
                Names of all reference sequences, should be unique.

            seqlengths:
                Lengths of all sequences in base pairs. This should contain
                non-negative values and have the same number of elements as
                ``seqnames``. Entries may also be None if no lengths are
                available for that sequence.

                Alternatively, a dictionary where keys are the sequence names
                and values are the lengths. If a name is missing from this
                dictionary, the length of the sequence is set to None.

                Alternatively a single integer, if all sequences are of the
                same length.

                Alternatively None, if no length information is available
                for any sequence.

            is_circular:
                Whether each sequence is circular. This should have the same
                number of elements as ``seqnames``. Entries may also be None
                if no information is available for that sequence.

                Alternatively, a dictionary where keys are the sequence names
                and values are the circular flags. If a name is missing from
                this dictionary, the flag for the sequence is set to None.

                Alternatively a single boolean, if all sequences have the same
                circular flag.

                Alternatively None, if no flags are available for any sequence.

            genome:
                The genome build containing each reference sequence. This
                should have the same number of elements as ``seqnames``.
                Entries may also be None if no information is available.

                Alternatively, a dictionary where keys are the sequence names
                and values are the genomes. If a name is missing from this
                dictionary, the genome is set to None.

                Alternatively a single string, if all sequences are derived
                from the same genome.

                Alternatively None, if no genome information is available
                for any sequence.

            validate:
                Whether to validate the arguments, internal use only.
        """
        self._seqnames = list(seqnames)
        self._reverse_seqnames = None
        self._seqlengths = self._flatten_incoming(seqlengths, int)
        self._is_circular = self._flatten_incoming(is_circular, bool)
        self._genome = self._flatten_incoming(genome, str)

        if validate:
            _validate_seqnames(self._seqnames)
            num_seqs = len(self._seqnames)
            _validate_seqlengths(self._seqlengths, num_seqs)
            _validate_is_circular(self._is_circular, num_seqs)
            _validate_genome(self._genome, num_seqs)

    def _populate_reverse_seqnames_index(self):
        if self._reverse_seqnames is None:
            revmap = {}
            for i, n in enumerate(self):
                if n not in revmap:
                    revmap[n] = i
            self._reverse_seqnames = revmap

    def _wipe_reverse_seqnames_index(self):
        self._reverse_seqnames = None

    def _flatten_incoming(self, values, expected) -> List:
        if values is None or isinstance(values, expected):
            return [values] * len(self)

        if isinstance(values, dict):
            output = []
            for n in self._seqnames:
                if n in values:
                    output.append(values[n])
                else:
                    output.append(None)
            return output

        values = _sanitize_vec(values)

        if isinstance(values, list):
            return values

        return list(values)

    def _define_output(self, in_place: bool = False) -> "SeqInfo":
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
            A deep copy of the current ``SeqInfo``.
        """
        from copy import deepcopy

        _seqnames_copy = deepcopy(self._seqnames)
        _seqlengths_copy = deepcopy(self._seqlengths)
        _is_circular_copy = deepcopy(self._is_circular)
        _genome_copy = deepcopy(self._genome)

        current_class_const = type(self)
        return current_class_const(
            seqnames=_seqnames_copy,
            seqlenghts=_seqlengths_copy,
            is_circular=_is_circular_copy,
            genome=_genome_copy,
            validate=False,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``SeqInfo``.
        """
        current_class_const = type(self)
        return current_class_const(
            self._seqnames,
            self._seqlengths,
            self._is_circular,
            self._genome,
            validate=False,
        )

    def copy(self):
        """Alias for :py:meth:`~__copy__`."""
        return self.__copy__()

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation of this ``SeqInfo``.
        """
        output = "SeqInfo(number_of_seqnames=" + str(len(self))
        output += ", seqnames=" + ut.print_truncated_list(self._seqnames)
        output += ", seqlengths=" + repr(self._seqlengths)
        output += ", is_circular=" + ut.print_truncated_list(self._is_circular)
        output += ", genome=" + ut.print_truncated_list(self._genome)

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this ``SeqInfo``.
        """
        output = f"SeqInfo with {len(self)} sequence{'s' if len(self) != 1 else ''}\n"

        nr = len(self)
        added_table = False
        if nr:
            if nr <= 10:
                indices = range(nr)
                insert_ellipsis = False
            else:
                indices = [0, 1, 2, nr - 3, nr - 2, nr - 1]
                insert_ellipsis = True

            raw_floating = ut.create_floating_names(None, indices)
            if insert_ellipsis:
                raw_floating = raw_floating[:3] + [""] + raw_floating[3:]
            floating = ["", ""] + raw_floating

            columns = []

            header = ["seqnames", "<str>"]
            showed = [f"{self._seqnames[x]}" for x in indices]
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            header = ["seqlengths", f"<{ut.print_type(self._seqlengths)}>"]
            showed = [f"{self._seqlengths[x]}" for x in indices]
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            header = ["is_circular", f"<{ut.print_type(self._is_circular)}>"]
            showed = [f"{self._is_circular[x]}" for x in indices]
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            header = ["genome", f"<{ut.print_type(self._genome)}>"]
            showed = [f"{self._genome[x]}" for x in indices]
            if insert_ellipsis:
                showed = showed[:3] + ["..."] + showed[3:]
            columns.append(header + showed)

            output += ut.print_wrapped_table(columns, floating_names=floating)
            added_table = True

        footer = []
        if len(footer):
            if added_table:
                output += "\n------\n"
            output += "\n".join(footer)

        return output

    ##########################
    ######>> seqnames <<######
    ##########################

    def get_seqnames(self) -> List[str]:
        """
        Returns:
            List of all chromosome names.
        """
        return self._seqnames

    def set_seqnames(
        self, seqnames: Sequence[str], in_place: bool = False
    ) -> "SeqInfo":
        """
        Args:
            seqnames:
                List of sequence names, of length equal to the number of names
                in this ``SeqInfo`` object. All names should be unique strings.

            in_place:
                Whether to modify the ``SeqInfo`` object in place.

        Returns:
            A modified ``SeqInfo`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _validate_seqnames(list(seqnames))

        output = self._define_output(in_place)
        output._seqnames = list(seqnames)
        return output

    @property
    def seqnames(self) -> List[str]:
        warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
        return self.get_seqnames()

    @seqnames.setter
    def seqnames(self, seqnames: Sequence[str]):
        warn(
            "Setting property 'seqnames' is an in-place operation, use 'set_seqnames' instead",
            UserWarning,
        )

        self.set_seqnames(seqnames, in_place=True)

    ############################
    ######>> seqlengths <<######
    ############################

    def get_seqlengths(self) -> List[int]:
        """
        Returns:
            A list of integers is returned containing the lengths of all
            sequences, in the same order as the sequence names from
            :py:meth:`~get_seqnames`.
        """
        return self._seqlengths

    def set_seqlengths(
        self,
        seqlengths: Optional[Union[int, Sequence[int], Dict[str, int]]],
        in_place: bool = False,
    ) -> "SeqInfo":
        """
        Args:
            seqlengths:
                List of sequence lengths, of length equal to the number of
                names in this ``SeqInfo`` object. Values may be None or
                non-negative integers.

                Alternatively, a dictionary where keys are the sequence
                names and values are the lengths. Not all names need to be
                present in which case the length is assumed to be None.

            in_place:
                Whether to modify the ``SeqInfo`` object in place.

        Returns:
            A modified ``SeqInfo`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        _seqlengths = self._flatten_incoming(seqlengths, int)
        _validate_seqlengths(_seqlengths, len(self))

        output = self._define_output(in_place)
        output._seqlengths = _seqlengths
        return output

    @property
    def seqlengths(self) -> List[int]:
        warn(
            "'seqlengths' is deprecated, use 'get_seqlengths' instead",
            UserWarning,
        )
        return self.get_seqlengths()

    @seqlengths.setter
    def seqlengths(
        self, seqlengths: Optional[Union[int, Sequence[int], Dict[str, int]]]
    ):
        warn(
            "Setting property 'seqlengths' is an in-place operation, use 'set_seqlengths' instead",
            UserWarning,
        )

        self.set_seqlengths(seqlengths, in_place=True)

    #############################
    ######>> is-circular <<######
    #############################

    def get_is_circular(self) -> List[bool]:
        """
        Returns:
            A list of booleans is returned specifying whether each sequence
            (from :py:meth:`~get_seqnames`) is circular.
        """
        return self._is_circular

    def set_is_circular(
        self,
        is_circular: Optional[Union[bool, Sequence[bool], Dict[str, bool]]],
        in_place: bool = False,
    ) -> "SeqInfo":
        """
        Args:
            is_circular:
                List of circular flags, of length equal to the number of
                names in this ``SeqInfo`` object. Values may be None or
                booleans.

                Alternatively, a dictionary where keys are the sequence
                names and values are the flags. Not all names need to be
                present in which case the flag is assumed to be None.

            in_place:
                Whether to modify the ``SeqInfo`` object in place.

        Returns:
            A modified ``SeqInfo`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _is_circular = self._flatten_incoming(is_circular, bool)
        _validate_is_circular(_is_circular, len(self))

        output = self._define_output(in_place)
        output._is_circular = _is_circular
        return output

    @property
    def is_circular(self) -> List[bool]:
        warn(
            "'is_circular' is deprecated, use 'get_is_circular' instead",
            UserWarning,
        )
        return self.get_is_circular()

    @is_circular.setter
    def is_circular(
        self, is_circular: Optional[Union[bool, Sequence[bool], Dict[str, bool]]]
    ):
        warn(
            "Setting property 'is_circular' is an in-place operation, use 'set_is_circular' instead",
            UserWarning,
        )

        self.set_is_circular(is_circular, in_place=True)

    ########################
    ######>> genome <<######
    ########################

    def get_genome(self) -> List[str]:
        """
        Returns:
            A list of strings is returned containing the genome identity for
            all sequences in :py:meth:`~get_seqnames`.
        """
        return self._genome

    def set_genome(
        self,
        genome: Optional[Union[str, Sequence[str], Dict[str, str]]],
        in_place: bool = False,
    ) -> "SeqInfo":
        """
        Args:
            genome:
                List of genomes, of length equal to the number of names in this
                ``SeqInfo`` object. Values may be None or strings.


            in_place:
                Whether to modify the ``SeqInfo`` object in place.

        Returns:
            A modified ``SeqInfo`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        _genome = self._flatten_incoming(genome, str)
        _validate_genome(_genome, len(self))

        output = self._define_output(in_place)
        output._genome = _genome
        return output

    @property
    def genome(self) -> List[str]:
        warn("'genome' is deprecated, use 'get_genome' instead", UserWarning)
        return self.get_genome()

    @genome.setter
    def genome(self, genome: Optional[Union[bool, Sequence[bool], Dict[str, bool]]]):
        warn(
            "Setting property 'genome' is an in-place operation, use 'set_genome' instead",
            UserWarning,
        )

        self.set_genome(genome, in_place=True)

    ######################################
    ######>> length and iterators <<######
    ######################################

    def __len__(self) -> int:
        """
        Returns:
            Number of sequences in this object.
        """
        return len(self._seqnames)

    def __iter__(self) -> SeqInfoIterator:
        """Iterator over sequences."""
        return SeqInfoIterator(self)

    #########################
    ######>> Slicers <<######
    #########################

    def get_subset(self, subset: Union[str, int, bool, Sequence]) -> "SeqInfo":
        """Subset ``SeqInfo``, based on their indices or seqnames.

        Args:
            subset:
                Indices to be extracted. This may be an integer, boolean, string,
                or any sequence thereof, as supported by
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.
                Scalars are treated as length-1 sequences.

                Strings may only be used if :py:attr:``~seqnames`` are available (see
                :py:meth:`~get_seqnames`). The first occurrence of each string
                in the seqnames is used for extraction.

        Returns:
            A new ``SeqInfo`` object with the sequences of interest.
        """
        if len(self) == 0:
            return SeqInfo.empty()

        idx, _ = ut.normalize_subscript(subset, len(self), self._seqnames)

        current_class_const = type(self)
        return current_class_const(
            seqnames=ut.subset_sequence(self._seqnames, idx),
            seqlengths=ut.subset_sequence(self._seqlengths, idx),
            is_circular=ut.subset_sequence(self._is_circular, idx),
            genome=ut.subset_sequence(self._genome, idx),
        )

    def __getitem__(self, subset: Union[str, int, bool, Sequence]) -> "SeqInfo":
        """Alias to :py:attr:`~get_subset`."""
        return self.get_subset(subset)

    @classmethod
    def empty(cls):
        """Create an zero-length `SeqInfo` object.

        Returns:
            same type as caller, in this case a `SeqInfo`.
        """
        return SeqInfo([], [], [], [])


@ut.combine_sequences.register
def _combine_SeqInfo(*x: SeqInfo) -> SeqInfo:
    return merge_SeqInfo(x)


def merge_SeqInfo(objects: List[SeqInfo]) -> SeqInfo:
    """Merge multiple :py:class:`~SeqInfo` objects, taking the union of all reference sequences. If the same reference
    sequence is present with the same details across ``objects``, only a single instance is present in the final object;
    if details are contradictory, they are replaced with None.

    Args:
        objects: List of ``SeqInfo`` objects.

    Returns:
        A single merged ``SeqInfo`` object.
    """
    all_sequences = {}

    for obj in objects:
        for i, y in enumerate(obj._seqnames):
            curlen = obj._seqlengths[i]
            curcir = obj._is_circular[i]
            curgen = obj._genome[i]

            if y not in all_sequences:
                all_sequences[y] = [curlen, curcir, curgen]
            else:
                present = all_sequences[y]
                prelen, precir, pregen = present
                if prelen != curlen:
                    present[0] = None
                if precir != curcir:
                    present[1] = None
                if pregen != curgen:
                    present[2] = None

    out_names = []
    out_lengths = []
    out_circular = []
    out_genome = []
    for k, v in all_sequences.items():
        out_names.append(k)
        out_lengths.append(v[0])
        out_circular.append(v[1])
        out_genome.append(v[2])

    return SeqInfo(out_names, out_lengths, out_circular, out_genome, validate=False)
