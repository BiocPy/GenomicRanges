from typing import Dict, List, Optional, Union, Sequence
import biocutils as bu
from warnings import warn


__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class SeqInfo:
    """Information about the reference sequences, specifically the name and length of each sequence, whether it is a
    circular, and the identity of the genome from which it was derived."""

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Optional[Union[int, Sequence[int], Dict[str, int]]],
        is_circular: Optional[Union[bool, Sequence[bool], Dict[str, bool]]],
        genome: Optional[Union[str, Sequence[str], Dict[str, str]]],
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
        self._seqlengths = self._flatten_incoming(seqlengths, int)
        self._is_circular = self._flatten_incoming(is_circular, bool)
        self._genome = self._flatten_incoming(genome, str)

        if validate:
            self._validate_seqnames()
            self._validate_seqlengths()
            self._validate_is_circular()
            self._validate_genome()

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
        if isinstance(values, list):
            return values
        return list(values)

    def _validate_seqnames(self):
        if not bu.is_list_of_type(self._seqnames, str):
            raise ValueError("'seqnames' should be a list of strings")
        n = len(self._seqnames)
        if n != len(set(self._seqnames)):
            raise ValueError("'seqnames' should contain unique strings")

    def _validate_seqlengths(self):
        n = len(self._seqnames)
        if not bu.is_list_of_type(self._seqlengths, int, ignore_none=True):
            raise ValueError("'seqlengths' should be a list of integers")
        if n != len(self._seqlengths):
            raise ValueError("'seqnames' and 'seqlengths' should have the same length")
        for l in self._seqlengths:
            if l < 0:
                raise ValueError("all entries of 'seqlengths' should be non-negative")

    def _validate_is_circular(self):
        n = len(self._seqnames)
        if not bu.is_list_of_type(self._is_circular, bool, ignore_none=True):
            raise ValueError("'is_circular' should be a list of booleans")
        if n != len(self._is_circular):
            raise ValueError("'seqnames' and 'is_circular' should have the same length")

    def _validate_genome(self):
        n = len(self._seqnames)
        if not bu.is_list_of_type(self._genome, str, ignore_none=True):
            raise ValueError("'genome' should be a list of strings")
        if n != len(self._genome):
            raise ValueError("'seqnames' and 'genome' should have the same length")

    @property
    def seqnames(self) -> List[str]:
        warn("'seqnames' is deprecated, use 'get_seqnames' instead", DeprecationWarning)
        return self.get_seqnames()

    def get_seqnames(self) -> List[str]:
        """
        Returns:
            List of all chromosome names.
        """
        return self._seqnames

    @property
    def seqlengths(self) -> List[int]:
        warn(
            "'seqlengths' is deprecated, use 'get_seqlengths' instead",
            DeprecationWarning,
        )
        return self.get_seqlengths()

    def get_seqlengths(self) -> List[int]:
        """
        Returns:
            A list of integers is returned containing the lengths of all
            sequences, in the same order as the sequence names from
            :py:meth:`~get_seqnames`.
        """
        return self._seqlengths

    @property
    def is_circular(self) -> List[bool]:
        warn(
            "'is_circular' is deprecated, use 'get_is_circular' instead",
            DeprecationWarning,
        )
        return self.get_is_circular()

    def get_is_circular(self) -> List[bool]:
        """
        Returns:
            A list of booleans is returned specifying whether each sequence
            (from :py:meth:`~get_seqnames`) is circular.
        """
        return self._is_circular

    @property
    def genome(self) -> List[str]:
        warn("'genome' is deprecated, use 'get_genome' instead", DeprecationWarning)
        return self.get_genome()

    def get_genome(self) -> List[str]:
        """
        Returns:
            A list of strings is returned containing the genome identity for
            all sequences in :py:meth:`~get_seqnames`.
        """
        return self._genome

    def _setter_copy(self, in_place: bool = False) -> "SeqInfo":
        if in_place:
            return self
        else:
            return type(self)(
                self._seqnames,
                self._seqlengths,
                self._is_circular,
                self._genome,
                validate=False,
            )

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
        output = self._setter_copy(in_place)
        output._seqnames = list(seqnames)
        output._validate_seqnames()
        return output

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
        output = self._setter_copy(in_place)
        output._seqlengths = self._flatten_incoming(seqlengths, int)
        output._validate_seqlengths()
        return output

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
        output = self._setter_copy(in_place)
        output._is_circular = self._flatten_incoming(is_circular, bool)
        output._validate_is_circular()
        return output

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
        output = self._setter_copy(in_place)
        output._genome = self._flatten_incoming(genome, str)
        output._validate_genome()
        return output

    def __len__(self) -> int:
        """
        Returns:
            Number of sequences in this object.
        """
        return len(self._seqnames)


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
