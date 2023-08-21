from itertools import zip_longest
from typing import Any, List, MutableMapping, Optional, Sequence, Union

from biocframe import BiocFrame

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class SeqInfo(BiocFrame):
    """Class that stores information about gene model.

    Must contain column "seqnames".

    Args:
        data (MutableMapping[str, Union[List[Any], MutableMapping]]): info about each
            sequence or chromosome. must contain a column `seqnames`.
        metadata (MutableMapping, optional): Additional metadata. Defaults to None.

    Raises:
        ValueError: If ``data`` does not contain required attributes.
    """

    required_columns = ["seqnames"]
    can_contain = ["seqnames", "seqlengths", "is_circular", "genome"]

    def __init__(
        self,
        data: MutableMapping[str, Union[List[Any], MutableMapping]],
        number_of_rows: Optional[int] = None,
        row_names: Optional[Sequence[str]] = None,
        column_names: Optional[Sequence[str]] = None,
        metadata: Optional[MutableMapping] = None,
    ) -> None:
        """Initialize a SeqInfo object."""
        super().__init__(data, number_of_rows, row_names, column_names, metadata)

    def _validate(self):
        """Internal function to validate SeqInfo."""

        if "genome" in self._data:
            if self._metadata is None:
                self._metadata = {}

            self._metadata["genome"] = self._data["genome"]
            del self._data["genome"]

        super()._validate()
        self._validate_seqs()

    def _validate_seqs(self):
        """Internal function to validate sequence information.

        Raises:
            ValueError: If missing required columns.
        """
        missing = list(set(self.required_columns).difference(set(self.column_names)))

        if len(missing) > 0:
            raise ValueError(
                f"data must contain {self.required_columns}. missing {missing}"
            )

    @property
    def seqnames(self) -> Sequence[str]:
        """Get sequence or chromosome names.

        Returns:
            Sequence[str]: List of all chromosome names.
        """
        return self.column("seqnames")

    @property
    def seqlengths(self) -> Optional[MutableMapping[str, int]]:
        """Get sequence or chromosome names and their lengths.

        Returns:
            (MutableMapping[str, int], optional): A dictionary containing chromosome names
                with their lengths.
        """

        if "seqlengths" not in self.column_names:
            return None

        return dict(zip_longest(self.column("seqnames"), self.column("seqlengths")))

    @property
    def is_circular(self) -> Optional[MutableMapping[str, bool]]:
        """Whether the sequences are circular.

        Returns:
            (MutableMapping[str, bool], optional): A dictionary containing chromosome names
                and if they are circular.
        """

        if "is_circular" not in self.column_names:
            return None

        return dict(zip_longest(self.column("seqnames"), self.column("is_circular")))

    @property
    def genome(self) -> Optional[str]:
        """Get genome/species information, if available.

        Returns:
            (str, optional): The species name or genome.
        """

        if self.metadata and "genome" in self.metadata:
            return self.metadata["genome"]

        return None

    @genome.setter
    def genome(self, genome: Optional[str]):
        """Set genome/species information.

        Args:
            genome (str): Species name or genome.
        """

        if self.metadata is None:
            self.metadata = {}

        self.metadata["genome"] = genome
