from itertools import zip_longest
from typing import Any, Dict, List, Optional

from biocframe import BiocFrame
from biocframe.types import DataType

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class SeqInfo(BiocFrame):
    """Class that stores information about gene model.

    Must contain column "seqnames".

    Args:
        data (Dict[str, Union[List[Any], Dict]]): Info about each
            sequence or chromosome. must contain a column `seqnames`.
        metadata (Dict, optional): Additional metadata. Defaults to None.

    Raises:
        ValueError: If ``data`` does not contain required attributes.
    """

    required_columns = ["seqnames"]
    can_contain = ["seqnames", "seqlengths", "is_circular", "genome"]

    def __init__(
        self,
        data: DataType,
        number_of_rows: Optional[int] = None,
        row_names: Optional[List[str]] = None,
        column_names: Optional[List[str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Initialize a SeqInfo object."""
        if "genome" in data:
            if metadata is None:
                metadata = {}

            metadata["genome"] = data["genome"]

            try:
                del data["genome"]
            except Exception:
                data = {k: v for k, v in data.items() if k != "genome"}

        super().__init__(
            data, number_of_rows, row_names, column_names, metadata
        )

        self._validate_seqs()

    def _validate_seqs(self):
        """Internal function to validate sequence information.

        Raises:
            ValueError: If missing required columns.
        """
        missing = list(
            set(self.required_columns).difference(set(self.column_names))
        )

        if len(missing) > 0:
            raise ValueError(
                f"data must contain {self.required_columns}. missing {missing}"
            )

    @property
    def seqnames(self) -> List[str]:
        """Get sequence or chromosome names.

        Returns:
            List[str]: List of all chromosome names.
        """
        return self.column("seqnames")

    @property
    def seqlengths(self) -> Optional[Dict[str, int]]:
        """Get sequence or chromosome names and their lengths.

        Returns:
            (Dict[str, int], optional): A dictionary containing chromosome names
                with their lengths.
        """

        if "seqlengths" not in self.column_names:
            return None

        return dict(
            zip_longest(self.column("seqnames"), self.column("seqlengths"))
        )

    @property
    def is_circular(self) -> Optional[Dict[str, bool]]:
        """Whether the sequences are circular.

        Returns:
            (Dict[str, bool], optional): A dictionary containing chromosome names
                and if they are circular.
        """

        if "is_circular" not in self.column_names:
            return None

        return dict(
            zip_longest(self.column("seqnames"), self.column("is_circular"))
        )

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
