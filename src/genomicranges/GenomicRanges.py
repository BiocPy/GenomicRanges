from ast import Num, Slice
from typing import Union, List, Dict
import pandas as pd
import ncls
import logging
from .ucsc import access_gtf_ucsc
from .gtf import parse_gtf
from .utils import split_pandas_df

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges:
    """GenomicRanges class to represent genomic locations"""

    def __init__(
        self,
        index: Dict[str, Union[ncls.NCLS32, ncls.NCLS64]],
        indices: List[Num],
        ranges: pd.DataFrame,
        metadata: Union[pd.DataFrame, None] = None,
    ) -> None:
        """Initialize an instance of GenomicRanges

        Args:
            index (Dict[str, Union[ncls.NCLS32, ncls.NCLS64]]): Nested Containment List index
            indices (List[Num]): indices that map between `index`, `ranges` and `metadata`
            ranges (pd.DataFrame): intervals as pandas dataframe
            metadata (Union[pd.DataFrame, None], optional): additional metadata as pandas dataframe. Defaults to None.
        """
        self.index = index
        self.indices = indices
        self.ranges = ranges
        self.metadata = metadata
        self._iterIdx = 0

    def seqnames(self) -> pd.Series:
        """access sequences or chromosome names

        Returns:
            pd.Series: sequence information across all positions
        """
        logging.info(self.ranges["seqnames"].value_counts())
        return self.ranges["seqnames"]

    def ranges(self) -> pd.DataFrame:
        """return the genomic positions

        Returns:
            pd.DataFrame: Dataframe containing the genomic positions
        """
        return self.ranges

    def strand(self) -> pd.Series:
        """access strand information (if available)

        Returns:
            pd.Series: strand across all columns
        """
        logging.info(self.ranges["strand"].value_counts())
        return self.ranges["strand"]

    def granges(self) -> "GenomicRanges":
        """GenomicRanges object with only ranges (no metadata)

        Returns:
            GenomicRanges: Genomic Ranges with only ranges
        """
        return GenomicRanges(self.index, self.indices, self.ranges)

    def mcols(self) -> pd.DataFrame:
        """Access metadata information about positions

        Returns:
            pd.DataFrame: metadata across all genomic positions
        """
        return self.metadata

    def len(self) -> Num:
        """Number of intervals

        Returns:
            Num: number of genomic internals
        """
        return self.ranges.shape[0]

    def length(self) -> Num:
        """Number of intervals

        Returns:
            Num: number of genomic internals
        """
        return self.len()

    def __len__(self):
        """Number of intervals

        Returns:
            Num: number of genomic internals
        """
        return self.len()

    def __iter__(self) -> "GenomicRanges":
        """Iterator"""
        return self

    def __next__(self) -> tuple:
        """Get value from current iteration

        Raises:
            StopIteration: when Index is out of range
        """
        try:
            if self._iterIdx < self.ranges.shape[0]:
                r = self.ranges.iloc[self._iterIdx]
                m = (
                    self.metadata.iloc[self._iterIdx]
                    if self.metadata is not None
                    else None
                )
        except IndexError:
            self._iterIdx = 0
            raise StopIteration()

        self._iterIdx += 1
        if self._iterIdx > self.ranges.shape[0]:
            raise StopIteration()

        return (r, m)

    def __str__(self) -> str:
        """string representation

        Returns:
            str: what the instance holds
        """
        return f"<GenomicRanges> contains {len(self.index.keys())} chromosomes, {self.ranges.shape[0]} intervals"

    def __getitem__(self, args: Union[Slice, List[int]]) -> "GenomicRanges":
        """Slice the object

        Args:
            args (Union[Slice, List[int]]): A slice object or a list of indices to subset

        Raises:
            Exception: args must be either a slice or list of indices less than the shape
            IndexError: indices exceeds the dataset dimensions

        Returns:
            GenomicRanges: a subset GenomicRanges object
        """
        argmin = float("-inf")
        argmax = float("inf")
        if isinstance(args, slice):
            argmin = args.start
            argmax = args.stop
        elif isinstance(args, list):
            argmin = min(args)
            argmax = max(args)
        else:
            raise Exception(
                f"{args} must be either a slice or list of indices."
            )

        if argmin == float("-inf") or argmax == float("inf"):
            raise Exception(f"slice not supported")

        # check if all indices are less than (min and max)
        for i in [argmin, argmax]:
            if i > self.ranges.shape[0]:
                logging.error(f"Indices must be between {self.ranges.shape}")
                return IndexError(f"index: {i} out of range")

        arg_slices = self.indices[args]

        df = self.ranges.iloc[arg_slices, :]
        if self.metadata is not None:
            df = pd.concat(
                [
                    df,
                    self.metadata.iloc[
                        arg_slices,
                    ],
                ],
                ignore_index=False,
                axis=1,
            )

        (indexes, indices, ranges, metadata) = split_pandas_df(df)
        return GenomicRanges(indexes, indices, ranges, metadata)

    def nearest(
        self, x: "GenomicRanges", k: int = 1
    ) -> List[Union[int, List[int]]]:
        """Find nearest neighbors that overlap with the positions in `x`

        Args:
            x (GenomicRanges): input genomic positions to find
            k (int, optional): find k nearest neighbors ?. Defaults to 1.

        Returns:
            List[Union[int, List[int]]]: Possible hits
        """
        hits = []
        for gr in x:
            grhits = self.index[gr[0]["seqnames"]].find_overlap(
                gr[0]["starts"], gr[0]["ends"]
            )
            counter = 0
            tmp_hits = []
            for i in grhits:
                if counter < k:
                    tmp_hits.append(i[2])
                else:
                    break
            hits.append(tmp_hits if len(tmp_hits) > 0 else None)
        return hits

    def toDF(self) -> pd.DataFrame:
        """Export or Convert `GenomicRanges` to Pandas DataFrame

        Returns:
            pd.DataFrame: Pandas DataFrame
        """
        df = self.ranges
        if self.metadata is not None:
            df = pd.concat([df, self.metadata], ignore_index=False, axis=1)

        return df

    @staticmethod
    def fromPandas(data: pd.DataFrame) -> "GenomicRanges":
        """Convert a pandas Dataframe to GenomicRanges

        Args:
            data (pd.DataFrame): a Pandas DataFrame object containing genomic positions.
                Must contain `seqnames`, `starts` & `ends` columns.

        Raises:
            Exception: Validation Error

        Returns:
            GenomicRanges: An Object representing genomic positions
        """

        # validation:
        # if `seqnames`, `starts` and `ends` don't exist, abort!
        if not set(["seqnames", "starts", "ends"]).issubset(
            set(data.columns.tolist())
        ):
            logging.error(
                f"DataFrame does not contain columns: `seqnames`, `starts` and `ends`"
            )
            raise Exception(
                f"DataFrame does not contain columns: `seqnames`, `starts` and `ends`"
            )

        (indexes, indices, ranges, metadata) = split_pandas_df(data)
        return GenomicRanges(indexes, indices, ranges, metadata)

    @staticmethod
    def fromGTF(file: str) -> "GenomicRanges":
        """Load a GTF file as GenomicRanges

        Args:
            file (str): path to gtf file

        Returns:
            GenomicRanges:  An Object representing genomic positions
        """
        compressed = True if file.endswith("gz") else False
        data = parse_gtf(file, compressed=compressed)

        (indexes, indices, ranges, metadata) = split_pandas_df(data)

        return GenomicRanges(indexes, indices, ranges, metadata)

    @staticmethod
    def fromUCSC(genome: str, type: str = "refGene") -> "GenomicRanges":
        """Load a GTF file from UCSC as GenomicRanges

        Args:
            genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
            type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

        Returns:
            GenomicRanges:  An Object representing genomic positions
        """
        path = access_gtf_ucsc(genome, type=type)
        compressed = True
        data = parse_gtf(path, compressed=compressed)

        (indexes, indices, ranges, metadata) = split_pandas_df(data)

        return GenomicRanges(indexes, indices, ranges, metadata)
