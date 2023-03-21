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
    adjust_interval,
)

import genomicranges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges(BiocFrame):
    """Class `GenomicRanges` represents genomic regions and annotations.
    
    Methods are available to load genome annotations from UCSC or 
    from a pandas DataFrame.

    Note: Pandas DataFrame must contain columns seqnames, starts, ends and strand.
    If strand column is not provided, * is used as the default value for each 
    genomic interval.

    ***Intervals are inclusive on both ends.***

    Additionally, `GenomicRanges` can also contain `Sequence Information` as 
    part of its metadata. In most cases, it usually contains for each sequence name 
    (or chromosome) in the object, its length. Sequence Information can also 
    contain additional metadata about the `genome`, and if its circular 
    (`isCircular`) or not.

    All columns other than `seqnames`, `starts, `ends` and `strand` are considered 
    metadata columns and can be accessed by the function `mcols`.

    Checkout the `SeqInfo` class for more information.

    Note: The documentation for some of the methods are copied over from the R 
        `GenomicRanges` package.
    """

    required_columns = ["seqnames", "starts", "ends", "strand"]

    def __init__(
        self,
        data: Mapping[str, Union[List[Any], Mapping]],
        numberOfRows: Optional[int] = None,
        rowNames: Optional[Sequence[str]] = None,
        columnNames: Optional[Sequence[str]] = None,
        metadata: Optional[Mapping] = None,
    ) -> None:
        """Initialize a `GenomicRanges` object.

        Note: `data` must contain `seqnames`, `starts` and `ends` columns.
        If `strand` column is not provided, `*` is used as the default value
        for each genomic interval.

        Args:
            data (Mapping[str, Union[List[Any], Mapping]]): 
                columns as dictionary, Must contain `seqnames`, `starts` and 
                `ends` columns.
            numberOfRows (int, optional): Number of genomic 
                intervals (or rows). Defaults to None.
            rowNames (Sequence[str], optional): Row index. 
                Defaults to None.
            columnNames (Sequence[str], optional): column names, 
                automatically inferred from `data`. Defaults to None.
            metadata (Mapping, optional): metadata. Defaults to None.
        """
        super().__init__(data, numberOfRows, rowNames, columnNames, metadata)

    def _validate(self):
        """Internal function to validate `GenomicRanges`.
        """
        if "strand" not in self._data:
            self._data["strand"] = ["*"] * len(self._data["starts"])

            if self._columnNames is not None:
                self._columnNames.append("strand")

        super()._validate()
        self._validate_ranges()

    def _validate_ranges(self):
        """Internal function to validate all columns of `GenomicRanges`.

        Raises:
            ValueError: If missing required columns.
        """
        missing = list(set(self.required_columns).difference(set(self.columnNames)))

        if len(missing) > 0:
            raise ValueError(
                f"data must contain {self.required_columns}."
                f"missing {missing} column{'s' if len(missing) > 1 else ''}"
            )

    @property
    def seqnames(self) -> Sequence[str]:
        """Get sequence or chromosome names.

        Returns:
            Sequence[str]: list of all chromosome names.
        """
        return self.column("seqnames")

    @property
    def start(self) -> Sequence[int]:
        """Get sequence or chromosome start positions.

        Returns:
            Sequence[int]: list of all chromosome start positions.
        """
        return self.column("starts")

    @property
    def end(self) -> Sequence[int]:
        """Get sequence or chromosome end positions.

        Returns:
            Sequence[int]: list of all chromosome end positions.
        """
        return self.column("ends")

    def ranges(
        self, ignoreStrand: bool = False, returnType: Optional[Callable] = None
    ) -> Union[pd.DataFrame, MutableMapping, "GenomicRanges", Any]:
        """Get genomic positions.

        Args:
            ignoreStrand (bool): ignore strands? Defaults to False.
            returnType (Callable, optional): format to return genomic positions. 
                Defaults to a dictionary representation, supports `pd.DataFrame`
                or any callable representation that takes a dictionary as an input.

        Raises:
            ValueError: `returnType` is not supported.

        Returns:
            Union[pd.DataFrame, MutableMapping, "GenomicRanges", Any]: genomic regions.
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
    def strand(self) -> Sequence[str]:
        """Get strand information. If strand is originally not provided, 
        we use '*' as a default value for each interval.

        Returns:
            Sequence[str]: strand across all positions.
        """
        return self.column("strand")

    @property
    def width(self) -> Sequence[int]:
        """Get widths of each interval.

        Returns:
            Sequence[int]: width of each interval.
        """

        widths = []

        for _, row in self:
            widths.append(row["ends"] - row["starts"])

        return widths

    @property
    def seqInfo(self) -> Optional["SeqInfo"]:
        """Get the sequence information object (if available).

        Returns:
            Optional["SeqInfo"]: Sequence information.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"]

        return None

    @seqInfo.setter
    def seqInfo(self, seqInfo: "SeqInfo"):
        """Set sequence information.

        Raises:
            ValueError: if property `seqInfo` is not a `SeqInfo` class object.

        Args:
            seqinfo ("SeqInfo"): sequence information.
        """

        if not isinstance(seqInfo, SeqInfo):
            raise ValueError("seqInfo is not a `SeqInfo` class object")

        if self._metadata is None:
            self._metadata = {}

        self._metadata["seqInfo"] = seqInfo

    @property
    def seqlengths(self) -> Optional[MutableMapping[str, int]]:
        """Get length of each chromosome (if available).

        Returns:
            Optional[MutableMapping[str, int]]: Sequence lengths.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].seqlengths

        return None

    @property
    def score(self) -> Optional[Sequence[Union[int, float]]]:
        """Get score (if available) for each genomic interval.

        Returns:
             Optional[Sequence[Union[int, float]]]: score column.
        """

        if "score" in self._columnNames:
            return self.data["score"]

        return None

    @score.setter
    def score(self, score: Sequence[Union[int, float]]):
        """Set score.

        Args:
            score (Sequence[Union[int, float]]): score values to set.

        Raises:
            ValueError: if length of provided `score` does not 
                match the number of intervals.
            TypeError: if `score` is not a list.
        """

        if not isinstance(score, list):
            raise TypeError("`score` must be a list!")

        if len(score) != self._numberOfRows:
            raise ValueError(
                "provided incorrect number of `score` values"
                f"must be {self._numberOfRows}, but provided {len(score)}"
            )

        self["score"] = score
        return self

    @property
    def isCircular(self) -> Optional[MutableMapping[str, bool]]:
        """Are the sequences/chromosomes circular? (only if available).

        Returns:
            Optional[MutableMapping[str, bool]]: a dictionary with 
                keys as chromosome and boolean value if its circular or not.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].isCircular

        return None

    @property
    def genome(self) -> Optional[str]:
        """Get genome information (if available).

        Returns:
            Optional[str]: genome.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].genome

        return None

    def granges(self) -> "GenomicRanges":
        """Creates a new `GenomicRanges` object with only ranges 
        (`seqnames`, `starts, `ends` and `strand`).

        Returns:
            GenomicRanges: `GenomicRanges` with only ranges.
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
        """Get metadata across all genomic intervals. 
        
        
        All columns other than `seqnames`, `starts, `ends` and `strand` 
        are considered metadata for each interval.

        Args:
            returnType (Callable, optional): format to return metadata. 
                Defaults to dictionary representation, supports `pd.DataFrame`
                or any callable representation that takes a dictionary as an input.

        Raises:
            ValueError: if `returnType` is not supported.

        Returns:
            Union[pd.DataFrame, MutableMapping, Any]: metadata columns without 
            genomic positions.
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
        pattern = (
            f"Class GenomicRanges with {self.dims[0]} intervals and "
            f"{self.dims[1] - 4} metadata columns \n"
            f"  columnNames: {self.columnNames}"
        )
        return pattern

    def __getitem__(
        self, args: Union[Sequence[str], Tuple[Sequence, Optional[Sequence]]]
    ) -> "GenomicRanges":
        """Slice a `GenomicRanges` object.

        Args:
            args (Union[Sequence[str], Tuple[Sequence, Optional[Sequence]]]): 
                indices to slice.
                if args is
                - Sequence[str]: a list of names, we slice by column names.
                - Tuple[Sequence, Optional[Sequence]]]: a tuple with indices, 
                    we slice by indices along the row and column axes.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with the subset.
        """
        new_frame = super().__getitem__(args)
        return GenomicRanges(
            new_frame._data,
            new_frame._numberOfRows,
            new_frame._rowNames,
            new_frame._columnNames,
            new_frame._metadata,
        )

    # intra-range methods
    def flank(
        self,
        width: int,
        start: bool = True,
        both: bool = False,
        ignoreStrand: bool = False,
    ) -> "GenomicRanges":
        """Generates flanking ranges for each range in the `GenomicRanges` 
        object. The logic for this - 

        (from the R/GenomicRanges & IRanges packages)

        - If `start` is `True` for a given range, the flanking occurs at the start, 
        otherwise the end. 
        - The `widths` of the flanks are given by the `width` parameter. 
        The widths can be negative, in which case the flanking region is 
        reversed so that it represents a prefix or suffix of the range. 

        Example:
            gr.flank(3, True), where x indicates a range in gr and 
            - indicates the resulting flanking region:

                ---xxxxxxx

            If start were FALSE, the range in gr becomes

                xxxxxxx---
            
            For negative width, i.e. gr.flank(x, -3, FALSE), 
                where * indicates the overlap between x and the result:

                xxxx***

            If both is True, then, for all ranges in x, 
                the flanking regions are extended into 
                (or out of, if width is negative) the range, 
                so that the result straddles the given endpoint 
                and has twice the width given by width. 
        
            This is illustrated below for gr.flank(3, both=TRUE):

                ---***xxxx

        Args:
            width (int): width to flank by.
            start (bool, optional): only flank starts?. Defaults to True.
            both (bool, optional): both starts and ends?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with the flanked ranges.
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
            if both is True:
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
        """Resize ranges to the specified `width` where either the `start`, 
        `end`, or `center` is used as an anchor.

        Args:
            width (int): width to resize, cannot be negative!.
            fix (str, optional): fix positions by `start`, `end` or `center`. 
                Defaults to "start".
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Raises:
            ValueError: if parameter fix is neither `start`, `end` or `center`.
            ValueError: if width is negative

        Returns:
            GenomicRanges: a new `GenomicRanges` object with the resized ranges.
        """

        if width < 0:
            raise ValueError("width cannot be negative!")

        if fix not in ["start", "end", "center"]:
            raise ValueError(
                f"`fix` must be either 'start', 'end' or 'resize', provided {fix}"
            )

        new_starts = []
        new_ends = []

        for idx, row in self:
            ts = None
            te = None

            if ignoreStrand is True or row["strand"] != "-":
                if fix == "start":
                    ts = row["starts"]
                    te = row["starts"] + width - 1
                elif fix == "center":
                    tmid = math.ceil((row["starts"] + row["ends"]) / 2)
                    twidthby2 = (
                        math.floor(width / 2)
                        if row["strand"] == "+"
                        else math.ceil(width / 2)
                    )
                    ts = tmid - twidthby2
                    te = ts + width - 1
                else:
                    te = row["ends"]
                    ts = row["ends"] - width + 1
            elif row["strand"] == "-":
                if fix == "end":
                    ts = row["starts"]
                    te = row["starts"] + width - 1
                elif fix == "center":
                    tmid = math.ceil((row["starts"] + row["ends"]) / 2)
                    twidthby2 = math.ceil(width / 2)
                    ts = tmid - twidthby2
                    te = ts + width - 1
                else:
                    te = row["ends"]
                    ts = row["ends"] - width + 1
            else:
                raise ValueError(
                    "strand must be either +, - or *, contains "
                    f"{row['strand']} at index: {idx}"
                )

            new_starts.append(ts)
            new_ends.append(te)

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
        """Shift all intervals by parameter `shift`.

        Args:
            shift (int, optional): shift interval. Defaults to 0.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with the shifted ranges.
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
        """Extend intervals to promoter regions.

        Generates promoter ranges relative to the transcription start site (TSS), 
        where TSS is start(x). The promoter range is expanded around the TSS 
        according to the upstream and downstream arguments. upstream represents 
        the number of nucleotides in the 5' direction and downstream the number 
        in the 3' direction. The full range is defined as, (start(x) - upstream) 
        to (start(x) + downstream - 1).

        Args:
            upstream (int, optional): number of positions to extend in the 5' 
                direction. Defaults to 2000.
            downstream (int, optional): number of positions to extend in the 3' 
                direction. Defaults to 200.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with ranges replaced 
            by promoter regions.
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
        """Restrict ranges to a given start and end positions.

        Args:
            start (int, optional): start position. Defaults to None.
            end (int, optional): end position. Defaults to None.
            keepAllRanges (bool, optional): Keep intervals that do 
                not overlap with start and end?. Defaults to False.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with restricted ranges.
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
        """Trim sequences outside of bounds for non-circular chromosomes.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with trimmed ranges.
        """

        # let just show a warning, shouldn't be an error
        warning_msg = """
        Not enough information to trim,
        `GenomicRanges` object does not contain sequence information.
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
        
        Important: these parameters are relative shift in positions for each range.

        Args:
            start (int, optional): relative start position. Defaults to None.
            width (int, optional): relative end position. Defaults to None.
            end (int, optional): relative width of the interval. Defaults to None.

        Raises:
            ValueError: if `width` is provided, either `start` or `end` must be provided.
            ValueError: provide two of the three parameters - `start`, `end` and `width` 
                but not all.

        Returns:
            GenomicRanges:  a new `GenomicRanges` object with narrow ranges.
        """
        if start is not None and end is not None and width is not None:
            raise ValueError(
                "only provide two of the three parameters - start, "
                "end and width but not all"
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
        """Internal method to calculate gap widths.

        Args:
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            Sequence[int]: gap widths for each range.
        """
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
        """Reduce orders the ranges, then merges overlapping or adjacent ranges.
        
        Args:
            withRevMap (bool, optional): return map of indices back to 
                original object?. Defaults to False.
            minGapwidth (int, optional): Ranges separated by a gap of 
                at least `minGapwidth` positions are not merged. Defaults to 1.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with reduced intervals.
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

        return genomicranges.fromPandas(finale)

    def range(
        self, withRevMap: bool = False, ignoreStrand: bool = False
    ) -> "GenomicRanges":
        """Calculate ranges for each chromosome.
        (minimum of all starts, maximum of all ends) in the object.
        
        Technically its same as `reduce` with a ridiculously high `minGapwidth`.

        Args:
            withRevMap (bool, optional): return map of indices back to 
                original object?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with the new ranges.
        """
        return self.reduce(
            withRevMap=withRevMap, minGapwidth=100000000000, ignoreStrand=ignoreStrand
        )

    def computeSeqlengths(self) -> MutableMapping[str, int]:
        """Get seqlengths either from the `SeqInfo` object 
            or computes one from the object.

        Note: if computed, they ae specific to this `GenomicRanges` 
        and may not represent the seqlenths of the genome.

        Returns:
            MutableMapping[str, int]: a dict of chromosome names and lengths.
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
        """Identify gaps in genomic positions for each distinct 
        `seqname` (chromosome) and `strand` combination.

        Args:
            start (int, optional): restrict chromosome start position. Defaults to 1.
            end (Optional[MutableMapping[str, int]], optional): restrict end 
                position for each chromosome. Defaults to None.
                If None, it uses the `SeqInfo` object if available.

        Returns:
            Optional["GenomicRanges"]: A new `GenomicRanges` containing 
            gaps across chromosome and strand. If there are no gaps, returns None.
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
        return genomicranges.fromPandas(final_df)

    def disjoin(
        self, withRevMap: bool = False, ignoreStrand: bool = False
    ) -> Optional["GenomicRanges"]:
        """Calculate disjoint genomic positions for each distinct 
        `seqname` (chromosome) and `strand` combination.

        Args:
            withRevMap (bool, optional):return map of indices back to 
                original object? . Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            Optional["GenomicRanges"]: A new `GenomicRanges` containing 
            disjoint ranges across chromosome and strand.
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
        return genomicranges.fromPandas(final_df)

    def isDisjoint(self, ignoreStrand: bool = False) -> bool:
        """Are all ranges (for each chromosome, strand) disjoint?

        Args:
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            bool: True if the intervals are disjoint.
        """
        new_gr = self.disjoin(ignoreStrand=ignoreStrand)

        return new_gr.shape != self.shape

    def disjointBins(self, ignoreStrand: bool = False):
        raise NotImplementedError

    # set operations
    def union(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find union of genomic intervals with `other`.

        Args:
            other (GenomicRanges): the other `GenomicRanges` object.

        Raises:
            TypeError: if other is not of type `GenomicRanges`.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with all ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("other is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

    def intersect(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find intersection of genomic intervals with `other`.

        Args:
            other (GenomicRanges): the other `GenomicRanges` object.

        Raises:
            TypeError: if other is not of type `GenomicRanges`.

        Returns:
            Optional["GenomicRanges"]: a new `GenomicRanges` object 
            with intersection ranges. If there are no 
            overlapping intervals, returns None.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("other is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

    def setdiff(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find set difference of genomic intervals with `other`.

        Args:
            other (GenomicRanges): the other `GenomicRanges` object.

        Raises:
            TypeError: if other is not of type `GenomicRanges`.
            
        
        Returns:
            Optional["GenomicRanges"]: a new `GenomicRanges` object 
            with the diff ranges. If there are intervals, returns None.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("other is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

    def binnedAverage(
        self, scorename: str, bins: "GenomicRanges", outname: str
    ) -> "GenomicRanges":
        """Calculate average for a column across all bins, then
        set a column called `outname` with those values.

        Args:
            scorename (str): the column to compute averages on.
            bins (GenomicRanges): bins you want to use.
            outname (str): new column name to add to the object.

        Raises:
            ValueError: if scorename column does not exist.
            Exception: scorename is not all ints or floats.
            TypeError: if bins is not of type `GenomicRanges`.

        Returns:
            GenomicRanges: a new `GenomicRanges` similar to bins, 
            with a new column containing the averages.
        """

        if not isinstance(bins, GenomicRanges):
            raise TypeError("bins is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

    def subtract(
        self, x: "GenomicRanges", minoverlap: int = 1, ignoreStrand: bool = False
    ):
        raise NotImplementedError

    # integer range methods
    def coverage(
        self, shift: int = 0, width: Optional[int] = None, weight: int = 1
    ) -> MutableMapping[str, np.ndarray]:
        """Calculate coverage for each chromosome, For each position, 
        counts the number of ranges that cover it.

        Args:
            shift (int, optional): shift all genomic positions. Defaults to 0.
            width (int, optional): restrict the width of all 
                chromosomes. Defaults to None.
            weight (int, optional): weight to use. Defaults to 1.

        Returns:
            MutableMapping[str, np.ndarray]: coverage vector for each chromosome.
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
        """Find overlaps between subject (self) and a query `GenomicRanges`.

        Args:
            query (GenomicRanges): query `GenomicRanges`.
            queryType (str, optional): overlap query type, must be one of 
                - "any": any overlap is good
                - "start": overlap at the beginning of the intervals
                - "end": must overlap at the end of the intervals
                - "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. 
                Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            Optional["GenomicRanges"]: A new `GenomicRanges` object 
            with the same length as query, containing hits to overlapping indices.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("query is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

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
            query (GenomicRanges): query `GenomicRanges`.
            queryType (str, optional): overlap query type, must be one of 
                - "any": any overlap is good
                - "start": overlap at the beginning of the intervals
                - "end": must overlap at the end of the intervals
                - "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. 
                Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            List[int]: number of overlaps for each range in query.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("query is not a `GenomicRanges` object")

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
            query (GenomicRanges): query `GenomicRanges`.
            queryType (str, optional): overlap query type, must be one of 
                - "any": any overlap is good
                - "start": overlap at the beginning of the intervals
                - "end": must overlap at the end of the intervals
                - "within": Fully contain the query interval. 
                Defaults to any.
            maxGap (int, optional): maximum gap allowed in the overlap. 
                Defaults to -1 (no gap allowed).
            minOverlap (int, optional): minimum overlap with query. Defaults to 1. 
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            Optional["GenomicRanges"]: A new `GenomicRanges` object 
            containing only subsets for overlaps in query.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("query is not a `GenomicRanges` object")

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
        """Internal function to search self and a query `GenomicRanges` object.

        Raises:
            TypeError: if query is not of type `GenomicRanges`
            
        Returns:
            Optional["GenomicRanges"]: a new `GenomicRanges` object that has 
            the same length as query but contains `hits` to 
            `indices` and `distance`.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("query is not a `GenomicRanges` object")

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
        return genomicranges.fromPandas(final_df)

    def nearest(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions both upstream and downstream that overlap with 
        each range in `query`. 
            
        Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): query `GenomicRanges` to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            "GenomicRanges": a new `GenomicRanges` containing list of 
            possible `hits` indices for each range in `query`.
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand)

    def precede(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with 
        each range in `query`.

        Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): query `GenomicRanges` to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            "GenomicRanges": a new `GenomicRanges` containing list of possible 
            `hits` indices for each range in `query`
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand, stepend=0)

    def follow(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only upstream that overlap with 
        each range in `query`. 
            
        Adds a new column to query called `hits`.

        Args:
            query (GenomicRanges): query `GenomicRanges` to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            "GenomicRanges": a new `GenomicRanges` containing list of possible 
            `hits` indices for each range in `query`.
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand, stepstart=0)

    def distanceToNearest(
        self, query: "GenomicRanges", ignoreStrand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with 
        each range in `query`. 
            
        Adds a new column to query called `hits`.
        
        Note: Technically same as nearest since we also return `distances`.

        Args:
            query (GenomicRanges): query `GenomicRanges` to find nearest positions.
            ignoreStrand (bool, optional): ignore strand? Defaults to False.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            "GenomicRanges": a new `GenomicRanges` containing list of possible 
            `hits` indices for each range in `query`.
        """
        return self._generic_search(query=query, ignoreStrand=ignoreStrand)

    # compare and order methods
    def duplicated(self,) -> Sequence[bool]:
        """Element wise comparison to find duplicate ranges.

        Returns:
            Sequence[bool]: True if duplicated else False.
        """
        df = self._generic_pandas_ranges(sort=False)
        return df.duplicated().to_list()

    def match(self, query: "GenomicRanges") -> Sequence[Optional[int]]:
        """Element wise comparison to find exact match ranges.

        Args:
            query (GenomicRanges): Input `GenomicRanges` to search matches.

        Raises:
            TypeError: if query is not of type `GenomicRanges`.
            
        Returns:
            Sequence[Optional[int]]: index if able to match else None.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("query is not a `GenomicRanges` object")

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
        """Internal function to create a pandas `DataFrame` from ranges.

        Args:
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.
            sort (bool, optional): sort by region?. Defaults to False.

        Returns:
            pd.DataFrame: a pandas `DataFrame` object.
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
        """Internal function to compute order.

        Args:
            query (GenomicRanges): query `GenomicRanges` to search matches.

        Returns:
            Sequence[Optional[int]]: sorted index order.
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
        """Sort the `GenomicRanges` object.

        Args:
            decreasing (bool, optional): decreasing order?. Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            "GenomicRanges": a new sorted `GenomicRanges` object.
        """
        order = self._generic_order(ignoreStrand=ignoreStrand)

        if decreasing:
            order = order[::-1]

        new_order = order.to_list()
        return self[new_order, :]

    def rank(self) -> Sequence[int]:
        """Get rank of the `GenomicRanges` object. 
        for each ranges identifies its position is a sorted order.

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
        """Split each range into chunks by `n` (number of sub intervals) 
        or `width` (intervals with equal width).
        
        Note: Either `n` or `width` must be provided, but not both.

        Args:
            n (int, optional): number of intervals to split into. 
                Defaults to None.
            width (int, optional): width of each interval. Defaults to None.

        Raises:
            ValueError: if both `n` and `width` are provided.

        Returns: 
            "GenomicRanges": a new `GenomicRanges` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either `n` or `width` must be provided but not both")

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
        return genomicranges.fromPandas(final_df)

    def tile(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each interval by `n` (number of sub intervals) 
        or `width` (intervals with equal width).
        
        Note: Either `n` or `width` must be provided but not both.

        Args:
            n (int, optional): number of intervals to split into. Defaults to None.
            width (int, optional): width of each interval. Defaults to None.

        Raises:
            ValueError: if both `n` and `width` are provided.

        Returns: 
            "GenomicRanges": a new `GenomicRanges` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("either `n` or `width` must be provided but not both")

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
        return genomicranges.fromPandas(final_df)

    def slidingWindows(self, width: int, step: int = 1) -> "GenomicRanges":
        """Slide along each range by `width` (intervals with equal `width`) and `step`.

        Args:
            width (Optional[int], optional): width of each interval. Defaults to None.
            step (int, optional), step interval, Defaults to 1.

        Returns: 
            "GenomicRanges": a new `GenomicRanges` with the sliding ranges.
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
        return genomicranges.fromPandas(final_df)

    # misc methods
    def sample(self, k: int = 5) -> "GenomicRanges":
        """Randomly sample `k` intervals from the object.

        Args:
            k (int, optional): number of intervals to sample. Defaults to 5.

        Returns:
            GenomicRanges: a new `GenomicRanges` with randomly sampled ranges.
        """
        sample = random.sample(range(len(self)), k=k)
        return self[sample, :]

    def invertStrand(self) -> "GenomicRanges":
        """Invert strand information for each interval.

        Conversion map: + becomes -, - becomes +, * stays the same.

        Returns:
            GenomicRanges: new `GenomicRanges` object.
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
        """Concatenate multiple `GenomicRanges` objects.

        Args:
            granges ("GenomicRanges"): `GenomicRanges` objects to concatenate. 

        Raises:
            TypeError: if any of the provided objects are not "GenomicRanges".

        Returns:
            GenomicRanges: new concatenated `GenomicRanges` object.
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
