from typing import (
    Union,
    List,
    Mapping,
    Any,
    Optional,
    Sequence,
    Tuple,
    MutableMapping,
    Callable,
)
from collections import OrderedDict
import pandas as pd
import numpy as np
import math
import random
import logging

from biocframe import BiocFrame
from .SeqInfo import SeqInfo
from .utils import (
    calc_row_gapwidth,
    find_disjoin,
    find_gaps,
    find_union,
    find_intersect,
    find_diff,
    compute_mean,
    create_np_interval_vector,
    OVERLAP_QUERY_TYPES,
    find_overlaps,
    find_nearest,
    split_intervals,
    slide_intervals,
)

from .ucsc import access_gtf_ucsc
from .gtf import parse_gtf

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges(BiocFrame):
    """Class GenomicRanges to represent genomic regions and annotations"""

    required_columns = ["seqnames", "starts", "ends", "strand"]

    def __init__(
        self,
        data: Optional[Mapping[str, Union[List[Any], Mapping]]] = ...,
        numberOfRows: Optional[int] = None,
        rowNames: Optional[Sequence[str]] = None,
        columnNames: Optional[Sequence[str]] = None,
        metadata: Optional[Mapping] = None,
    ) -> None:
        super().__init__(data, numberOfRows, rowNames, columnNames, metadata)

        # self._setIndex()

    # @staticmethod
    # def buildIndex(
    #     seqnames: Union[Sequence[str], pd.Series],
    #     starts: Union[Sequence[int], pd.Series],
    #     ends: Union[Sequence[int], pd.Series],
    # ) -> Mapping[str, Union[ncls.NCLS32, ncls.NCLS64]]:
    #     """Given genomic positions, builds an NCLS index

    #     Args:
    #         seqnames (Union[Sequence[str], pd.Series]): sequence or chromosome names
    #         starts (Union[Sequence[int], pd.Series]): genomic start interval
    #         ends (Union[Sequence[int], pd.Series]): genomic end interval

    #     Returns:
    #         Tuple[Mapping[str, Union[ncls.NCLS32, ncls.NCLS64]], pd.DataFrame]: a tuple containing
    #             an NCLS for each chromsome and a pandas Dataframe of all ranges.
    #     """
    #     ranges = pd.DataFrame({"seqnames": seqnames, "starts": starts, "ends": ends})
    #     ranges["_index"] = range(0, ranges.shape[0])
    #     groups = ranges.groupby("seqnames")

    #     # generate NCLS indexes for each seqname
    #     indexes = {}
    #     for group, rows in groups:
    #         indexes[group] = ncls.NCLS(
    #             rows.starts.astype(int).values,
    #             rows.ends.astype(int).values,
    #             rows._index.astype(int).values,
    #         )

    #     return indexes

    # def _setIndex(self):
    #     """Internal function to set or update NCLS index
    #     """
    #     self._index = GenomicRanges.buildIndex(
    #         self.column("seqnames"), self.column("starts"), self.column("ends")
    #     )

    def _validate(self):
        """internal function to validate GenomicRanges
        """

        if "strand" not in self._data:
            self._data["strand"] = ["*"] * len(self._data["starts"])

            if self._columnNames is not None:
                self._columnNames.append("strand")

        super()._validate()
        self._validate_ranges()

    def _validate_ranges(self):
        """Internal function to validate genomic ranges

        Raises:
            ValueError: if missing required columns
        """
        missing = list(set(self.required_columns).difference(set(self.columnNames)))

        if len(missing) > 0:
            raise ValueError(
                f"data must contain {self.required_columns}. missing {missing}"
            )

    @property
    def seqnames(self) -> Sequence[str]:
        """Get sequence or chromosome names

        Returns:
            Sequence[str]: list of all chromosome names
        """
        return self.column("seqnames")

    @property
    def start(self) -> Sequence[int]:
        """Get sequence or chromosome start positions

        Returns:
            Sequence[int]: list of all chromosome start positions
        """
        return self.column("starts")

    @property
    def end(self) -> Sequence[int]:
        """Get sequence or chromosome end positions

        Returns:
            Sequence[int]: list of all chromosome end positions
        """
        return self.column("ends")

    def ranges(
        self, ignoreStrand: bool = False, returnType: Optional[Callable] = None
    ) -> Union[pd.DataFrame, MutableMapping, "GenomicRanges"]:
        """Get the genomic positions

        Args:
            ignoreStrand (bool): ignore strands? Defaults to False.
            returnType (Callable, optional): format to return genomic positions. 
                Defaults to dictionary structure. currently supports `pd.DataFrame`.

        Raises:
            ValueError: `returnType` is not supported

        Returns:
            Union[pd.DataFrame, MutableMapping, "GenomicRanges"]: genomic regions
        """

        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
        }

        if ignoreStrand:
            obj["strand"] = ["*"] * len(obj["seqnames"])

        if returnType is None:
            return obj
        else:
            try:
                return returnType(obj)
            except Exception as e:
                raise ValueError(f"{returnType} not supported, {str(e)}")

    @property
    def strand(self) -> Optional[Sequence[str]]:
        """Get strand information (if available)

        Returns:
            Optional[Sequence[str]]: strand across all positions or None
        """
        return self.column("strand")

    @property
    def width(self) -> Sequence[int]:
        """Get widths of each interval

        Returns:
            Sequence[int]: width of each interval
        """

        widths = []

        for _, row in self:
            widths.append(row["ends"] - row["starts"])

        return widths

    @property
    def seqInfo(self) -> Optional["SeqInfo"]:
        """Get the sequence information object (if available)

        Returns:
            SeqInfo: Sequence information
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"]

        return None

    @seqInfo.setter
    def seqInfo(self, seqInfo: "SeqInfo"):
        """Set sequence information

        Raises:
            ValueError: seqInfo is a `SeqInfo` class object

        Args:
            seqinfo ("SeqInfo"): sequence information
        """

        if not isinstance(seqInfo, SeqInfo):
            raise ValueError("seqInfo is not a `SeqInfo` class object")

        if self._metadata is None:
            self._metadata = {}

        self._metadata["seqInfo"] = seqInfo

    @property
    def seqlengths(self) -> Optional[MutableMapping[str, int]]:
        """Get length of each chromosome (if available)

        Returns:
            Optional[MutableMapping[str, int]]: Sequence lengths
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].seqlengths

        return None

    @property
    def score(self) -> Optional[Sequence[Union[int, float]]]:
        """Get the score (if available).

        Returns:
             Optional[Sequence[Union[int, float]]]: score
        """

        if "score" in self._columnNames:
            return self.data["score"]

        return None

    @score.setter
    def score(self, score: Sequence[Union[int, float]]):
        """Get the score (if available).

        Args:
            score (Sequence[Union[int, float]]): score values to set.

        Returns:
             Optional[Sequence[Union[int, float]]]: score
        """

        self["score"] = score
        return self

    @property
    def isCircular(self) -> Optional[MutableMapping[str, bool]]:
        """are the sequences circular?

        Returns:
            Optional[MutableMapping[str, bool]]: a dict for each chromosome 
                and if its circular or not.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].isCircular

        return None

    @property
    def genome(self) -> Optional[str]:
        """Genome information.

        Returns:
            Optional[str]: the genome
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].genome

        return None

    def granges(self) -> "GenomicRanges":
        """New GenomicRanges object with only ranges (`seqnames`, `starts and `ends`)

        Returns:
            GenomicRanges: Genomic Ranges with only ranges
        """
        return GenomicRanges(
            {
                "seqnames": self.column("seqnames"),
                "starts": self.column("starts"),
                "ends": self.column("ends"),
                "strand": self.column("strand"),
            }
        )

    def mcols(
        self, returnType: Optional[Callable] = None
    ) -> Union[pd.DataFrame, MutableMapping]:
        """Access metadata about positions

        Args:
            returnType (Callable, optional): format to return metadata. 
                Defaults to dictionary structure. currently supports `pd.DataFrame`.

        Raises:
            ValueError: `returnType` is not supported

        Returns:
            Union[pd.DataFrame, MutableMapping]: metadata columns without positions
        """

        new_data = OrderedDict()
        for k in self.columnNames:
            if k not in self.required_columns:
                new_data[k] = self.column(k)

        if returnType is None:
            return new_data
        else:
            try:
                return returnType(new_data)
            except Exception as e:
                raise ValueError(f"{returnType} not supported, {str(e)}")

    def __str__(self) -> str:
        pattern = """
        Class GenomicRanges with {} intervals and {} metadata columns
          columnNames: {}
        """
        return pattern.format(self.dims[0], self.dims[1] - 3, self.columnNames)

    def __getitem__(
        self, args: Union[Sequence[str], Tuple[Sequence, Optional[Sequence]]]
    ) -> "GenomicRanges":
        """Slice a GenomicRanges object

        Args:
            args (Union[Sequence[str], Tuple[Sequence, Optional[Sequence]]]): indices to slice.
                Sequence[str]: Slice by column names
                Tuple[Sequence, Optional[Sequence]]]: slice by indices along the row and column axes.

        Returns:
            GenomicRanges: returns a new GenomicRanges object with the subset
        """
        new_frame = super().__getitem__(args)
        return GenomicRanges(
            new_frame._data,
            new_frame._numberOfRows,
            new_frame._rowNames,
            new_frame._columnNames,
            new_frame._metadata,
        )

    def _adjustInterval(
        self, row: MutableMapping[str, Any], shiftStart: int = 0, shiftEnd: int = 0,
    ) -> Tuple[int, int]:
        """internal function to shift genomic intervals

        Args:
            row (MutableMapping[str, Any]): a row from GenomicRanges
            shiftStart (int, optional): number of positions to shift start by. Defaults to 0.
            shiftEnd (int, optional): number of positions to shift end by. Defaults to 0.

        Returns:
            adjusted genomic intervals
        """
        return (row["starts"] + shiftStart, row["ends"] + shiftEnd)

    # intra-range methods
    def flank(
        self,
        width: int,
        start: bool = True,
        both: bool = False,
        ignoreStrand: bool = False,
    ) -> "GenomicRanges":
        """Recover regions flanking the set of ranges represented by the GenomicRanges object

        Args:
            width (int): width to flank by
            start (bool, optional): only flank starts?. Defaults to True.
            both (bool, optional): both starts and ends?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new GenomicRanges object with the flanked intervals
        """
        new_starts = []
        new_ends = []

        all_starts = self.column("starts")
        all_ends = self.column("ends")
        all_strands = self.column("strand")

        # figure out which position to pin, start or end?
        start_flags = [start] * len(all_strands)
        if not ignoreStrand:
            start_flags = [
                start != (all_strands[i] == "-") for i in range(len(all_strands))
            ]

        # if both is true, then depending on start_flag, we extend it out
        # I couldn't understand the scenario's with witdh <=0,
        # so refer to the R implementation here
        for idx in range(len(start_flags)):
            sf = start_flags[idx]
            tstart = 0
            if both:
                tstart = (
                    all_starts[idx] - abs(width)
                    if sf
                    else all_ends[idx] - abs(width) + 1
                )
            else:
                if width >= 0:
                    tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1
                else:
                    tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1

            new_starts.append(tstart)
            new_ends.append(tstart + (width * (2 if both else 1) - 1))

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def resize(
        self, width: int, fix: str = "start", ignoreStrand: bool = False,
    ) -> "GenomicRanges":
        """Resize intervals based on strand in the `GenomicRanges` object

        Args:
            width (int): width to resize
            fix (str, optional): fix positions by `start` or `end`. Defaults to "start".
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Raises:
            ValueError: parameter fix is neither `start` nor `end`

        Returns:
            GenomicRanges: a new GenomicRanges object with the resized intervals
        """

        if fix not in ["start", "end"]:
            raise ValueError(f"`fix` must be either 'start' or 'end', provided {fix}")

        new_starts = []
        new_ends = []

        all_starts = self.column("starts")
        all_ends = self.column("ends")
        all_strands = self.column("strand")

        # figure out which position to pin, start or end?
        start_flags = [fix == "start"] * len(all_strands)
        if not ignoreStrand:
            start_flags = [
                start_flags[i] != (all_strands[i] == "-")
                for i in range(len(all_strands))
            ]

        for idx in range(len(start_flags)):
            sf = start_flags[idx]
            tstart = all_starts[idx] if sf else all_ends[idx] - width + 1

            new_starts.append(tstart)
            new_ends.append(tstart + width - 1)

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def shift(self, shift: int = 0) -> "GenomicRanges":
        """Shift all intervals

        Args:
            shift (int, optional): shift interval. Defaults to 0.

        Returns:
            GenomicRanges: a new GenomicRanges object with the shifted intervals
        """
        if shift == 0:
            return self

        all_starts = [x + shift for x in self.column("starts")]
        all_ends = [x + shift for x in self.column("ends")]

        new_data = self._data.copy()
        new_data["starts"] = all_starts
        new_data["ends"] = all_ends

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def promoters(self, upstream: int = 2000, downstream: int = 200) -> "GenomicRanges":
        """Extend intervals to promoter regions

        Args:
            upstream (int, optional): number of positions to extend in the 5' direction . Defaults to 2000.
            downstream (int, optional): number of positions to extend in the 3' direction. Defaults to 200.

        Returns:
            GenomicRanges: a new GenomicRanges object with the extended intervals for promoter regions
        """
        all_starts = self.column("starts")
        all_ends = self.column("ends")
        all_strands = self.column("strand")

        start_flags = [all_strands[i] != "-" for i in range(len(all_strands))]

        new_starts = [
            (
                all_starts[idx] - upstream
                if start_flags[idx]
                else all_ends[idx] - downstream + 1
            )
            for idx in range(len(start_flags))
        ]
        new_ends = [
            (
                all_starts[idx] + downstream - 1
                if start_flags[idx]
                else all_ends[idx] + upstream
            )
            for idx in range(len(start_flags))
        ]

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def restrict(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        keepAllRanges: bool = False,
    ) -> "GenomicRanges":
        """Restrict intervals to a given start and end positions across all chormosomes

        Args:
            start (Optional[int], optional): start position. Defaults to None.
            end (Optional[int], optional): end position. Defaults to None.
            keepAllRanges (bool, optional): Keep intervals that do not overlap with start and end?. Defaults to False.

        Returns:
            GenomicRanges: a new GenomicRanges object with restricted intervals
        """
        all_starts = self.column("starts")
        all_ends = self.column("ends")

        new_data = self._data.copy()

        new_starts = all_starts.copy()
        if start:
            new_starts = [(start if x < start else x) for x in all_starts]

        new_ends = all_ends.copy()
        if end:
            new_ends = [(end if x > end else x) for x in all_ends]

        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        new_rowNames = self._rowNames.copy()

        if not keepAllRanges:
            keepers = []
            # find intervals indices that overlap with the given start and end
            for idx in range(len(all_starts)):
                keep = False
                if start and (all_starts[idx] <= start <= all_ends[idx]):
                    keep = True

                if end and (all_starts[idx] <= end <= all_ends[idx]):
                    keep = True

                if keep:
                    keepers.append(idx)

            for col in self._columnNames:
                new_data[col] = [new_data[col][x] for x in keepers]

            if self._rowNames:
                new_rowNames = [new_rowNames[x] for x in keepers]

        return GenomicRanges(
            new_data,
            columnNames=self._columnNames,
            rowNames=new_rowNames,
            metadata=self._metadata,
        )

    def trim(self) -> "GenomicRanges":
        """Trim sequences outside of bounds for non-circular chromosomes

        Returns:
            GenomicRanges: a new GenomicRanges object with trimmed positions
        """

        # let just show a warning, shouldn't be an error
        warning_msg = """
        Not enough information to trim,
        GenomicRanges object does not contain sequence information.
        """
        if not self.seqInfo:
            logging.warning(warning_msg)
            return self

        if self._metadata is None:
            logging.warning(warning_msg)
            return self

        if self._metadata["seqInfo"] is None:
            logging.warning(warning_msg)
            return self

        seqinfos = self.seqInfo
        seqlengths = seqinfos.seqlengths
        isCircular = seqinfos.isCircular

        if seqlengths is None:
            logging.warning(warning_msg)
            return self

        if isCircular is None:
            logging.warning("considering all sequences as non-circular...")

        all_chrs = self.column("seqnames")
        all_ends = self.column("ends")

        new_data = self._data.copy()
        new_rowNames = self._rowNames.copy()

        keepers = []
        for idx in range(len(all_chrs)):
            keep = True
            t_chr = all_chrs[idx]
            if (
                isCircular is not None
                and isCircular[t_chr] is False
                and all_ends[idx] > seqlengths[t_chr]
            ):
                keep = False

            if keep:
                keepers.append(idx)

        for col in self._columnNames:
            new_data[col] = [new_data[col][x] for x in keepers]

        if self._rowNames:
            new_rowNames = [new_rowNames[x] for x in keepers]

        return GenomicRanges(
            new_data,
            columnNames=self._columnNames,
            rowNames=new_rowNames,
            metadata=self._metadata,
        )

    # TODO: needs checks when relative - {start, width and end} do not agree
    def narrow(
        self,
        start: Optional[int] = None,
        width: Optional[int] = None,
        end: Optional[int] = None,
    ) -> "GenomicRanges":
        """Narrow genomic positions by provided start, width and end parameters. 
            Important: these parameters are relative shift in positions for each interval.

        Args:
            start (Optional[int], optional): relative start position. Defaults to None.
            width (Optional[int], optional): relative end position. Defaults to None.
            end (Optional[int], optional): relative width of the interval. Defaults to None.

        Raises:
            ValueError: when parameters were set incorrectly

        Returns:
            GenomicRanges:  a new GenomicRanges object with narrow positions
        """
        if start is not None and end is not None and width is not None:
            raise ValueError(
                "only provide two of the three parameters - start, end and width but not all"
            )

        if width is not None:
            if start is None and end is None:
                raise ValueError(
                    "if width is provided, either start or end must be provided"
                )

        new_starts = []
        new_ends = []

        all_starts = self.column("starts")
        all_ends = self.column("ends")

        if start:
            new_starts = [x + start - 1 for x in all_starts]
        else:
            new_starts = all_starts.copy()

        if end:
            new_ends = [x + end - 1 for x in all_starts]
        else:
            new_ends = all_ends

        if width:
            if start is None:
                new_starts = [x - width + 1 for x in new_ends]
            elif end is None:
                new_ends = [x + width - 1 for x in new_starts]
            else:
                raise ValueError(
                    "if width is provided, either start or end must be provided"
                )

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def _calcGapwidths(self, ignoreStrand: bool = False) -> Sequence[int]:
        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        df_gr = pd.DataFrame(obj)
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])

        if ignoreStrand:
            obj["strand"] = ["*"] * len(self.column("seqnames"))

        gapwidths = (
            df_gr["index"]
            .rolling(2)
            .apply(
                lambda x: calc_row_gapwidth(
                    df_gr.loc[x.index[0], :], df_gr.loc[x.index[1], :]
                )
            )
        )

        return gapwidths

    # inter range methods

    # TODO: implement dropEmptyRanges
    # TODO: this is a very ineffecient implementation, someone can do a better job later.
    def reduce(
        self,
        withRevMap: bool = False,
        minGapwidth: int = 1,
        ignoreStrand: bool = False,
    ) -> "GenomicRanges":
        """Reduce orders the ranges, then merges overlapping or adjacent intervals.
        
        Args:
            withRevMap (bool, optional): return map of indices back to original object?. Defaults to False.
            minGapwidth (int, optional): Ranges separated by a gap of at least minGapwidth positions are not merged. Defaults to 1.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new GenomicRanges object with reduced intervals
        """

        df_gr = self._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)

        df_gr["gapwidths"] = self._calcGapwidths(ignoreStrand=ignoreStrand)
        df_gr["gapwidth_flag"] = [
            (True if (pd.isna(x) or (x == 0)) else x < minGapwidth)
            for x in df_gr["gapwidths"]
        ]

        gaps_to_merge = df_gr[df_gr["gapwidth_flag"] == True]

        gaps_merged = gaps_to_merge.groupby(
            ["seqnames", "strand", "gapwidth_flag"], sort=False
        ).agg(starts=("starts", min), ends=("ends", max), revmap=("index", list))

        gaps_merged = gaps_merged.reset_index()

        gaps_not_merged = df_gr[df_gr["gapwidth_flag"] == False]
        gaps_not_merged["revmap"] = gaps_not_merged["index"].apply(lambda x: [x])
        gaps_not_merged = gaps_not_merged[gaps_merged.columns]

        finale = pd.concat([gaps_merged, gaps_not_merged])

        columns_to_keep = ["seqnames", "strand", "starts", "ends"]

        if withRevMap:
            columns_to_keep.append("revmap")

        finale = finale[columns_to_keep].sort_values(
            ["seqnames", "strand", "starts", "ends"]
        )

        return GenomicRanges.fromPandas(finale)

    def range(
        self, withRevMap: bool = False, ignoreStrand: bool = False
    ) -> "GenomicRanges":
        """Calculate ranges for each chromosome 
            (minimum of all starts, maximum of all ends) in the object.

            Technically its same as `reduce` with a ridiculously high minGapwidth.

        Args:
            withRevMap (bool, optional): return map of indices back to original object?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new GenomicRanges object with the new intervals
        """
        return self.reduce(
            withRevMap=withRevMap, minGapwidth=100000000000, ignoreStrand=ignoreStrand
        )

    def _computeSeqLengths(self) -> MutableMapping[str, int]:
        """Internal function to access seqlengths either from the SeqInfo object 
            or computes one by looking at the genomic positions.

        Returns:
            MutableMapping[str, int]: a dict of chromosome names and lengths
        """
        seqlengths = self.seqlengths

        if seqlengths is None:
            seqlengths = {}
            ranges = self.range()
            for _, row in ranges:
                seqlengths[row["seqnames"]] = row["ends"]

        return seqlengths

    def gaps(
        self, start: int = 1, end: Optional[MutableMapping[str, int]] = None
    ) -> Optional["GenomicRanges"]:
        """Identify gaps in genomic positions for each distinct chromosome and strand combination.

        Args:
            start (int, optional): restrict chromosome start position. Defaults to 1.
            end (Optional[MutableMapping[str, int]], optional): restrict end position for each chromosome. Defaults to None.
                If None, this tried to use the SeqInfo object if available.

        Returns:
            Optional["GenomicRanges"]: A new GenomicRanges containing gap regions across chromosome and strand
        """

        # seqlengths = self._computeSeqLengths()
        seqlengths = self.seqlengths

        if end is None:
            end = seqlengths

        df_gr = self._generic_pandas_ranges(sort=True)
        groups = df_gr.groupby(["seqnames", "strand"])

        gap_intervals = []
        for name, group in groups:
            group["iindex"] = range(len(group))

            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            end_limit = None
            if end is not None:
                if name[0] in end and end[name[0]] is not None:
                    end_limit = end[name[0]]

            tgaps = find_gaps(all_intvals, start_limit=start, end_limit=end_limit)

            for td in tgaps:
                tmp_start = td[0]
                if tmp_start < start:
                    tmp_start = start

                tmp_end = td[1]
                if end_limit is not None and tmp_end > end_limit:
                    tmp_end = end_limit

                td_res = (name[0], name[1], tmp_start, tmp_end)

                gap_intervals.append(td_res)

        if end is not None:
            groups = df_gr.groupby(["seqnames"])
            for name, group in groups:
                strands = group["strand"].unique()
                missing_strands = list(set(["*", "+", "-"]).difference(set(strands)))
                for ms in missing_strands:
                    gap_intervals.append((name, ms, start, end[name]))

        if len(gap_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(gap_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def disjoin(
        self, withRevMap: bool = False, ignoreStrand: bool = False
    ) -> Optional["GenomicRanges"]:
        """Calculate disjoint genomic positions for each distinct chromosome and strand combination.

        Args:
            withRevMap (bool, optional):return map of indices back to original object? . Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            Optional["GenomicRanges"]: A new GenomicRanges containing disjoint ranges across chromosome and strand
        """

        df_gr = self._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)
        groups = df_gr.groupby(["seqnames", "strand"])

        disjoin_intervals = []
        for name, group in groups:
            group["iindex"] = range(len(group))

            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            tdisjoin = find_disjoin(all_intvals, withRevMap=withRevMap)

            for td in tdisjoin:
                td_res = (name[0], name[1], td[0], td[1])
                if withRevMap:
                    td_res.append(td[2])
                disjoin_intervals.append(td_res)

        if len(disjoin_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        if withRevMap:
            columns.append("revmap")

        final_df = pd.DataFrame.from_records(disjoin_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def isDisjoint(self, ignoreStrand: bool = False) -> bool:
        """Are all ranges (for each chromosome, strand) disjoint?

        Args:
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            bool: True if the intervals are disjoint
        """
        new_gr = self.disjoin(ignoreStrand=ignoreStrand)

        return new_gr.shape != self.shape

    # def disjointBins(self, ignoreStrand: bool = False):
    #     raise NotImplementedError

    # set operations
    def union(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find union of genomic intervals with `other`

        Args:
            other (GenomicRanges): the other GenomicRanges object

        Returns:
            GenomicRanges: a new GenomicRanges object with 
                the intervals
        """
        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        df_gr = pd.concat([a_df_gr, b_df_gr])
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])
        groups = df_gr.groupby(["seqnames", "strand"])

        union_intervals = []
        for name, group in groups:
            group["iindex"] = range(len(group))

            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            tunion = find_union(all_intvals)

            for td in tunion:
                td_res = (name[0], name[1], td[0], td[1])
                union_intervals.append(td_res)

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(union_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def intersect(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find intersection of genomic intervals with `other`

        Args:
            other (GenomicRanges): the other GenomicRanges object

        Returns:
            Optional["GenomicRanges"]: a new GenomicRanges object with 
                the intervals
        """
        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        df_gr = pd.concat([a_df_gr, b_df_gr])
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])
        groups = df_gr.groupby(["seqnames", "strand"])

        intersect_intervals = []
        for name, group in groups:
            group["iindex"] = range(len(group))

            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            tintersect = find_intersect(all_intvals)

            for td in tintersect:
                td_res = (name[0], name[1], td[0], td[1])
                intersect_intervals.append(td_res)

        if len(intersect_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(intersect_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def setdiff(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find set difference of genomic intervals with `other`

        Args:
            other (GenomicRanges): the other GenomicRanges object

        Returns:
            Optional["GenomicRanges"]: a new GenomicRanges object with 
                the diff intervals
        """
        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        only_seqnames = pd.concat(
            [a_df_gr[["seqnames", "strand"]], b_df_gr[["seqnames", "strand"]]]
        )
        only_seqnames["seqstrand"] = only_seqnames["seqnames"] + only_seqnames["strand"]

        unique_seqs = list(only_seqnames["seqstrand"].unique())

        diff_ints = []
        for name in unique_seqs:
            ustrand = name[-1]
            chrom = name[:-1]

            a_set = a_df_gr[
                (a_df_gr["seqnames"] == chrom) & (a_df_gr["strand"] == ustrand)
            ]

            if len(a_set) == 0:
                continue

            a_intvals = [
                (x[0], x[1])
                for x in zip(a_set["starts"].to_list(), a_set["ends"].to_list())
            ]

            b_set = b_df_gr[
                (b_df_gr["seqnames"] == chrom) & (b_df_gr["strand"] == ustrand)
            ]
            b_intvals = [
                (x[0], x[1])
                for x in zip(b_set["starts"].to_list(), b_set["ends"].to_list())
            ]

            tdiff = find_diff(a_intvals, b_intvals)
            for td in tdiff:
                td_res = (chrom, ustrand, td[0], td[1])
                diff_ints.append(td_res)

        if len(diff_ints) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(diff_ints, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def binnedAverage(
        self, scorename: str, bins: "GenomicRanges", outname: str
    ) -> "GenomicRanges":
        """Calculate average for a score column across all bins, then
            set a column called outname with those values

        Args:
            scorename (str): the column to compute averages on
            bins (GenomicRanges): bins you want to use
            outname (str): new column name to add to bins

        Raises:
            ValueError: if scorename column does not exist
            Exception: scorename is not all ints or floats

        Returns:
            GenomicRanges: a new GenomicRanges with a column containing
                the averages
        """

        if scorename not in self.columnNames:
            raise ValueError(f"{scorename} is not a valid column name")

        values = self.column(scorename)

        if not all((isinstance(x, int) or isinstance(x, float)) for x in values):
            raise Exception(
                f"{scorename} is not a valid column, its neither ints not floats"
            )

        df_gr = self._generic_pandas_ranges(sort=True)
        df_gr["values"] = values

        tgt_gr = bins._generic_pandas_ranges(ignoreStrand=True, sort=True)
        tgt_groups = tgt_gr.groupby(["seqnames"])

        result = []
        cache_intvals = {}

        for name, group in tgt_groups:
            src_intervals = df_gr[df_gr["seqnames"] == name]

            if len(src_intervals) == 0:
                for _, g in group.iterrows():
                    result.append((name, g["starts"], g["ends"], "*", None))
                continue

            if name not in cache_intvals:
                all_intvals = [
                    (x[0], x[1])
                    for x in zip(
                        src_intervals["starts"].to_list(),
                        src_intervals["ends"].to_list(),
                    )
                ]

                _, np_sum = compute_mean(
                    intervals=all_intvals, values=src_intervals["values"].to_list()
                )

                cache_intvals[name] = np_sum

            np_sum = cache_intvals[name]

            for _, tint in group.iterrows():
                vec = np_sum[tint["starts"] - 1 : tint["ends"]]
                vec_mean = np.sum(vec) / np.count_nonzero(vec)

                result.append((name, tint["starts"], tint["ends"], "*", vec_mean))

        columns = ["seqnames", "starts", "ends", "strand", outname]
        final_df = pd.DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    # def subtract(
    #     self, x: "GenomicRanges", minoverlap: int = 1, ignoreStrand: bool = False
    # ):
    #     pass

    # integer range methods
    def coverage(
        self, shift: int = 0, width: Optional[int] = None, weight: int = 1
    ) -> MutableMapping[str, np.ndarray]:
        """Calculate coverage for each chromosome

        Args:
            shift (int, optional): shift all genomic positions. Defaults to 0.
            width (Optional[int], optional): restrict the width of all chromosomes. Defaults to None.
            weight (int, optional): weight to use. Defaults to 1.

        Returns:
            MutableMapping[str, np.ndarray]: _description_
        """
        df_gr = self._generic_pandas_ranges(sort=True)
        groups = df_gr.groupby(["seqnames"])

        shift_arr = None
        if shift > 0:
            shift_arr = np.zeros(shift)

        result = {}
        for name, group in groups:
            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            cov, _ = create_np_interval_vector(intervals=all_intvals, withRevMap=False)

            if shift > 0:
                cov = np.concatenate((shift_arr, cov))

            if weight > 0:
                cov = cov * weight

            if width is not None:
                cov = cov[:width]

            result[name] = cov

        return result

    # search based methods
    def findOverlaps(
        self,
        query: "GenomicRanges",
        queryType: str = "any",
        maxGap: int = -1,
        minOverlap: int = 1,
        ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Find overlaps between subject (self) and a query genomic ranges.

        Args:
            query (GenomicRanges): query genomic ranges.
            queryType (str, optional): overlap query type, must be one of 
                    "any": any overlap is good
                    "start": overlap at the beginning of the intervals
                    "end": must overlap at the end of the intervals
                    "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            Optional["GenomicRanges"]: A GenomicRanges object with the same length as query, 
                containing hits to overlapping indices.
        """

        if queryType not in OVERLAP_QUERY_TYPES:
            raise ValueError(f"{queryType} must be one of {OVERLAP_QUERY_TYPES}")

        df_gr = self._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)
        tgt_gr = query._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)

        tgt_groups = tgt_gr.groupby(["seqnames", "strand"])

        result = []
        for name, group in tgt_groups:
            chrom = name[0]
            ustrand = name[1]
            src_intervals = df_gr[
                (df_gr["seqnames"] == chrom) & (df_gr["strand"] == ustrand)
            ]

            src_intvals_map = src_intervals["index"].to_list()

            if len(src_intervals) == 0:
                for _, g in group.iterrows():
                    result.append((chrom, ustrand, g["starts"], g["ends"], []))
                continue

            subject_intvals = [
                (x[0], x[1])
                for x in zip(
                    src_intervals["starts"].to_list(), src_intervals["ends"].to_list()
                )
            ]

            query_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            thits = find_overlaps(
                subject_intvals,
                query_intvals,
                maxGap=maxGap,
                minOverlap=minOverlap,
                queryType=queryType,
            )

            for th in thits:
                tindices = [src_intvals_map[i - 1] for i in th[2]]
                result.append((chrom, ustrand, th[0][0], th[0][1], tindices))

        if len(result) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends", "hits"]
        final_df = pd.DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def countOverlaps(
        self,
        query: "GenomicRanges",
        queryType: str = "any",
        maxGap: int = -1,
        minOverlap: int = 1,
        ignoreStrand: bool = False,
    ) -> List[int]:
        """Count overlaps between `subject` (self) and a `query` genomic ranges.

        Args:
            query (GenomicRanges): query genomic ranges.
            queryType (str, optional): overlap query type, must be one of 
                    "any": any overlap is good
                    "start": overlap at the beginning of the intervals
                    "end": must overlap at the end of the intervals
                    "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            List[int]: number of overlaps for each interval in query
        """
        result = self.findOverlaps(
            query=query,
            queryType=queryType,
            maxGap=maxGap,
            minOverlap=minOverlap,
            ignoreStrand=ignoreStrand,
        )

        if result is None:
            return None

        hits = result.column("hits")
        return [len(ht) for ht in hits]

    def subsetByOverlaps(
        self,
        query: "GenomicRanges",
        queryType: str = "any",
        maxGap: int = -1,
        minOverlap: int = 1,
        ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Subset `subject` (self) with overlaps in `query` genomic ranges.

        Args:
            query (GenomicRanges): query genomic ranges.
            queryType (str, optional): overlap query type, must be one of 
                    "any": any overlap is good
                    "start": overlap at the beginning of the intervals
                    "end": must overlap at the end of the intervals
                    "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            Optional["GenomicRanges"]: A GenomicRanges object containings only subsets 
                for overlaps in query.
        """
        result = self.findOverlaps(
            query=query,
            queryType=queryType,
            maxGap=maxGap,
            minOverlap=minOverlap,
            ignoreStrand=ignoreStrand,
        )

        if result is None:
            return None

        hits = result.column("hits")
        hit_counts = [len(ht) for ht in hits]
        indices = [idx for idx in range(len(hit_counts)) if hit_counts[idx] > 0]

        return query[indices, :]

    def _generic_search(
        self,
        query: "GenomicRanges",
        ignoreStrand: bool = False,
        stepstart: int = 3,
        stepend: int = 3,
    ) -> Optional["GenomicRanges"]:
        """Internal function to search self and a query genomic ranges object

        Returns:
            Optional["GenomicRanges"]: a new GenomicRanges object that has the same length as query
                but contains `hits` to indices and `distance`.
        """
        subject_gr = self._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)
        query_gr = query._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)

        query_groups = query_gr.groupby(["seqnames", "strand"])

        result = []
        for name, group in query_groups:
            chrom = name[0]
            ustrand = name[1]
            src_intervals = subject_gr[
                (subject_gr["seqnames"] == chrom) & (subject_gr["strand"] == ustrand)
            ]

            src_intvals_map = src_intervals["index"].to_list()

            if len(src_intervals) == 0:
                for _, g in group.iterrows():
                    result.append((chrom, ustrand, g["starts"], g["ends"], [], None))
                continue

            subject_intvals = [
                (x[0], x[1])
                for x in zip(
                    src_intervals["starts"].to_list(), src_intervals["ends"].to_list()
                )
            ]

            query_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            thits = find_nearest(
                subject_intvals, query_intvals, stepstart=stepstart, stepend=stepend
            )
            for th in thits:
                tindices = [src_intvals_map[i - 1] for i in th[2]]
                result.append((chrom, ustrand, th[0][0], th[0][1], tindices, th[3]))

        if len(result) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends", "hits", "distance"]
        final_df = pd.DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def nearest(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions both upstream and downstream that overlap with the 
            each genomics interval in `query`. Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): input GenomicRanges to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Returns:
            "GenomicRanges": List of possible hit indices for each interval in `query`
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand)

    def precede(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with the 
            each genomics interval in `query`. Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): input GenomicRanges to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Returns:
            "GenomicRanges": List of possible hit indices for each interval in `query`
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand, stepend=0)

    def follow(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only upstream that overlap with the 
            each genomics interval in `query`. Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): input GenomicRanges to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Returns:
            "GenomicRanges": List of possible hit indices for each interval in `query`
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand, stepstart=0)

    def distanceToNearest(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with the 
            each genomics interval in `query`. Adds a new column to query called `hits`.
            Technically same as nearest since we also return `distances`.

        Args:
            query (GenomicRanges): input GenomicRanges to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Returns:
            "GenomicRanges": List of possible hit indices for each interval in `query`
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand)

    # compare and order methods
    def duplicated(self,) -> Sequence[bool]:
        """Element wise comparison to find duplicate intervals

        Returns:
            Sequence[bool]: True if duplicated else False
        """
        df = self._generic_pandas_ranges(sort=False)
        return df.duplicated().to_list()

    def match(self, query: "GenomicRanges") -> Sequence[Optional[int]]:
        """Element wise comparison to exact match intervals.

        Args:
            query (GenomicRanges): input GenomicRanges to search matches.

        Returns:
            Sequence[Optional[int]]: index if able to match else None
        """
        df = self._generic_pandas_ranges(sort=False)

        hits = []
        for _, row in query:
            sliced = df[
                (df["seqnames"] == row["seqnames"])
                & (df["strand"] == row["strand"])
                & (df["starts"] == row["starts"])
                & (df["ends"] == row["ends"])
            ]

            if len(sliced) == 0:
                hits.append(None)
            else:
                hits.append(list(sliced["index"].unique()))
        return hits

    def _generic_pandas_ranges(self, ignoreStrand=False, sort=False) -> pd.DataFrame:
        """Internal function to create a pandas dataframe from ranges

        Args:
            ignoreStrand (bool, optional): Ignore strand?. Defaults to False.
            sort (bool, optional): sort by region?. Defaults to False.

        Returns:
            pd.DataFrame: a Pandas DataFrame object
        """
        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        if ignoreStrand:
            obj["strand"] = ["*"] * len(self.column("seqnames"))

        df = pd.DataFrame(obj)

        if sort:
            df = df.sort_values(["seqnames", "strand", "starts", "ends"])

        return df

    def _generic_order(self, ignoreStrand=False) -> Sequence[int]:
        """Internal function that provides order

        Args:
            query (GenomicRanges): input GenomicRanges to search matches.

        Returns:
            Sequence[Optional[int]]: sorted index order
        """
        sorted = self._generic_pandas_ranges(ignoreStrand=ignoreStrand, sort=True)
        new_order = sorted["index"]

        return new_order

    def isUnsorted(self, ignoreStrand=False) -> bool:
        """Are the genomic positions unsorted?

        Args:
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            bool: True if unsorted else False.
        """
        order = self._generic_order(ignoreStrand=ignoreStrand)
        diff = order.diff()

        if any(diff != 1):
            return True

        return False

    def order(self, decreasing=False) -> Sequence[int]:
        """Get the order of indices for sorting.

        Args:
            decreasing (bool, optional): descending order?. Defaults to False.

        Returns:
            Sequence[int]: order of indices.
        """
        order = self._generic_order()

        if decreasing:
            order = order[::-1]
        return order.to_list()

    def sort(
        self, decreasing: bool = False, ignoreStrand: bool = False
    ) -> "GenomicRanges":
        """Sort the GenomicRanges object.

        Args:
            decreasing (bool, optional): decreasing order?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            "GenomicRanges": a new sorted GenomicRanges object.
        """
        order = self._generic_order(ignoreStrand=ignoreStrand)

        if decreasing:
            order = order[::-1]

        new_order = order.to_list()
        return self[new_order, :]

    def rank(self) -> Sequence[int]:
        """Get rank of the GenomicRanges object. 
            for each interval identifies its position is a sorted order

        Returns:
            Sequence[int]: list of indices identifying rank.
        """
        order = self._generic_order().to_list()
        rank = [order.index(x) for x in range(len(order))]
        return rank

    # windowing functions
    def tileByRange(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each range into chunks by 
                n (number of sub intervals) or width (intervals with equal width).
            Either n or width must be provided but not both.

        Args:
            n (Optional[int], optional): number of intervals to split into. Defaults to None.
            width (Optional[int], optional): width of each interval. Defaults to None.

        Raises:
            ValueError: when both n and width are both provided.

        Returns: 
            "GenomicRanges": a new GenomicRanges with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either n or width must be provided but not both")

        ranges = self.range()

        all_intervals = []
        for _, val in ranges:
            twidth = None
            if n is not None:
                twidth = int((val["ends"] - val["starts"]) / n)
            elif width is not None:
                twidth = width

            all_intervals.extend(
                split_intervals(
                    val["seqnames"], val["strand"], val["starts"], val["ends"], twidth
                )
            )

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def tile(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each interval by 
                n (number of sub intervals) or width (intervals with equal width).
            Either n or width must be provided but not both.

        Args:
            n (Optional[int], optional): number of intervals to split into. Defaults to None.
            width (Optional[int], optional): width of each interval. Defaults to None.

        Raises:
            ValueError: when both n and width are both provided.

        Returns: 
            "GenomicRanges": a new GenomicRanges with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either n or width must be provided but not both")

        all_intervals = []
        for _, val in self:
            twidth = None
            if n is not None:
                twidth = math.ceil((val["ends"] - val["starts"] + 1) / (n))
            elif width is not None:
                twidth = width

            all_intervals.extend(
                split_intervals(
                    val["seqnames"], val["strand"], val["starts"], val["ends"], twidth
                )
            )

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def slidingWindows(self, width: int, step: int = 1) -> "GenomicRanges":
        """slide along each interval by width (intervals with equal width) and step.

        Args:
            width (Optional[int], optional): width of each interval. Defaults to None.
            step (int, optional), step interval, Defaults to 1.

        Returns: 
            "GenomicRanges": a new GenomicRanges with the sliding ranges.
        """
        all_intervals = []
        for _, val in self:
            all_intervals.extend(
                slide_intervals(
                    val["seqnames"],
                    val["strand"],
                    val["starts"],
                    val["ends"],
                    width=width,
                    step=step,
                )
            )

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    # misc methods
    def sample(self, k: int = 5) -> "GenomicRanges":
        """Randomly sample `k` intervals from the Object.

        Args:
            k (int, optional): number of intervals to sample. Defaults to 5.

        Returns:
            GenomicRanges: a new GenomicRanges object.
        """
        sample = random.sample(range(len(self)), k=k)
        return self[sample, :]

    def invertStrand(self) -> "GenomicRanges":
        """Invert strand information for each interval.

        Returns:
            GenomicRanges: new GenomicRanges object.
        """
        convertor = {"+": "-", "-": "+", "*": "*"}
        inverts = [convertor[idx] for idx in self.column("strand")]

        new_data = self._data.copy()
        new_data["strand"] = inverts

        return GenomicRanges(
            new_data,
            numberOfRows=self._numberOfRows,
            rowNames=self._rowNames,
            columnNames=self._columnNames,
            metadata=self._metadata,
        )

    def concat(self, *granges: "GenomicRanges") -> "GenomicRanges":
        """Concatenate GenomicRanges objects.

        Args:
            granges ("GenomicRanges"): "GenomicRanges" to concatenate. 

        Raises:
            TypeError: if any of the provided objects are not "GenomicRanges".

        Returns:
            GenomicRanges: new concatenated "GenomicRanges" object.
        """
        all_granges = [isinstance(gr, GenomicRanges) for gr in granges]

        if not all(all_granges):
            raise TypeError("all provided objects are not GenomicRanges objects")

        all_columns = [gr.columnNames for gr in granges]
        all_columns.append(self.columnNames)
        all_unique_columns = list(
            set([item for sublist in all_columns for item in sublist])
        )

        new_data = OrderedDict()
        for col in all_unique_columns:
            if col not in new_data:
                new_data[col] = []

            for gr in granges:
                if col in gr.columnNames:
                    new_data[col].extend(gr.column(col))
                else:
                    new_data[col].extend([None] * len(gr))

            if col in self.columnNames:
                new_data[col].extend(self.column(col))
            else:
                new_data[col].extend([None] * len(self))

        return GenomicRanges(new_data, columnNames=all_unique_columns)

    @staticmethod
    def tileGenome(
        seqlengths: MutableMapping, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Create new genomic regions by partitioning a specified genome.

        Args:
            seqlengths (MutableMapping): sequence lengths of each chromosome.
            n (Optional[int], optional): number of intervals to split into. Defaults to None.
            width (Optional[int], optional): width of each interval. Defaults to None.

        Raises:
            ValueError: either n or width must be provided but not both

        Returns:
            GenomicRanges: a new GenomicRanges with the tile regions.
        """
        if n is not None and width is not None:
            raise ValueError("either n or width must be provided but not both")

        all_intervals = []
        for key, val in seqlengths.items():
            twidth = None
            if n is not None:
                twidth = math.ceil(val / n)
            elif width is not None:
                twidth = width

            all_intervals.extend(split_intervals(key, "*", 1, val, twidth))

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    @staticmethod
    def fromPandas(data: pd.DataFrame) -> "GenomicRanges":
        """Convert a pandas Dataframe to GenomicRanges

        Args:
            data (pd.DataFrame): a Pandas DataFrame object containing genomic positions.
                Must contain `seqnames`, `starts` & `ends` columns.

        Returns:
            GenomicRanges: An Object representing genomic positions
        """

        obj = OrderedDict()

        for col in data.columns:
            obj[col] = data[col].to_list()

        rindex = None
        if data.index is not None:
            rindex = data.index.to_list()

        return GenomicRanges(obj, rowNames=rindex, columnNames=data.columns.to_list())

    @staticmethod
    def fromGTF(file: str) -> "GenomicRanges":
        """Load a genome annotation from GTF file as `GenomicRanges`

        Args:
            file (str): path to gtf file

        Returns:
            GenomicRanges:  An Object representing genomic positions
        """
        compressed = True if file.endswith("gz") else False
        data = parse_gtf(file, compressed=compressed)

        return GenomicRanges.fromPandas(data)

    @staticmethod
    def fromUCSC(genome: str, type: str = "refGene") -> "GenomicRanges":
        """Load a genome annotation from UCSC as `GenomicRanges`

        Args:
            genome (str): genome shortcode; e.g. hg19, hg38, mm10 etc
            type (str): One of refGene, ensGene, knownGene or ncbiRefSeq

        Returns:
            GenomicRanges:  Gene model represented as GenomicRanges
        """
        path = access_gtf_ucsc(genome, type=type)
        compressed = True
        data = parse_gtf(path, compressed=compressed)

        return GenomicRanges.fromPandas(data)

