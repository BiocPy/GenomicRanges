from typing import (
    Union,
    List,
    Any,
    Optional,
    Sequence,
    MutableMapping,
)
from itertools import zip_longest
from biocframe import BiocFrame

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class SeqInfo(BiocFrame):
    required_columns = ["seqnames"]
    can_contain = ["seqnames", "seqlengths", "isCircular", "genome"]

    def __init__(
        self,
        data: MutableMapping[str, Union[List[Any], MutableMapping]],
        numberOfRows: Optional[int] = None,
        rowNames: Optional[Sequence[str]] = None,
        columnNames: Optional[Sequence[str]] = None,
        metadata: Optional[MutableMapping] = None,
    ) -> None:
        super().__init__(data, numberOfRows, rowNames, columnNames, metadata)

    def _validate(self):
        """Internal function to validate SeqInfo.
        """

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
            ValueError: if missing required columns.
        """
        missing = list(set(self.required_columns).difference(set(self.columnNames)))

        if len(missing) > 0:
            raise ValueError(
                f"data must contain {self.required_columns}. missing {missing}"
            )

    @property
    def seqnames(self) -> Sequence[str]:
        """Get sequence or chromosome names.

        Returns:
            Sequence[str]: list of all chromosome names.
        """
        return self.column("seqnames")

    @property
    def seqlengths(self) -> Optional[MutableMapping[str, int]]:
        """Get sequence or chromosome names and their lengths.

        Returns:
            MutableMapping[str, int]: dict containing chromosome names 
                with their lengths.
        """

        if "seqlengths" not in self._columnNames:
            return None

        return dict(zip_longest(self.column("seqnames"), self.column("seqlengths")))

    @property
    def isCircular(self) -> Optional[MutableMapping[str, bool]]:
        """are the sequences Circular?

        Returns:
            MutableMapping[str, bool]: dict containing chromosome names 
                and if they are circular.
        """

        if "isCircular" not in self._columnNames:
            return None

        return dict(zip_longest(self.column("seqnames"), self.column("isCircular")))

    @property
    def genome(self) -> Optional[str]:
        """Get genome/species information.

        Returns:
            str: get species name or genome.
        """

        if self._metadata and "genome" in self._metadata:
            return self._metadata["genome"]

        return None

    @genome.setter
    def genome(self, genome: str):
        """Set genome/species information

        Args:
            genome (str): species name or genome
        """

        if self._metadata is None:
            self._metadata = {}

        self._metadata["genome"] = genome
