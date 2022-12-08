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
import ncls

from biocframe import BiocFrame
from .SeqInfo import SeqInfo
from .utils import (
    calc_row_gapwidth,
    find_disjoin,
    find_gaps,
    find_union,
    find_intersect,
    find_diff,
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

        self._setIndex()

    @staticmethod
    def buildIndex(
        seqnames: Union[Sequence[str], pd.Series],
        starts: Union[Sequence[int], pd.Series],
        ends: Union[Sequence[int], pd.Series],
    ) -> Mapping[str, Union[ncls.NCLS32, ncls.NCLS64]]:
        """Given genomic positions, builds an NCLS index

        Args:
            seqnames (Union[Sequence[str], pd.Series]): sequence or chromosome names
            starts (Union[Sequence[int], pd.Series]): genomic start interval
            ends (Union[Sequence[int], pd.Series]): genomic end interval

        Returns:
            Tuple[Mapping[str, Union[ncls.NCLS32, ncls.NCLS64]], pd.DataFrame]: a tuple containing
                an NCLS for each chromsome and a pandas Dataframe of all ranges.
        """
        ranges = pd.DataFrame({"seqnames": seqnames, "starts": starts, "ends": ends})
        ranges["_index"] = range(0, ranges.shape[0])
        groups = ranges.groupby("seqnames")

        # generate NCLS indexes for each seqname
        indexes = {}
        for group, rows in groups:
            indexes[group] = ncls.NCLS(
                rows.starts.astype(int).values,
                rows.ends.astype(int).values,
                rows._index.astype(int).values,
            )

        return indexes

    def _setIndex(self):
        """Internal function to set or update NCLS index
        """
        self._index = GenomicRanges.buildIndex(
            self.column("seqnames"), self.column("starts"), self.column("ends")
        )

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

        # let just print a warning, shouldn't be an error
        warning_msg = """
        Not enough information to trim,
        GenomicRanges object does not contain sequence information.
        """
        if not self.seqInfo:
            print(warning_msg)
            return self

        if self._metadata is None:
            print(warning_msg)
            return self

        if self._metadata["seqInfo"] is None:
            print(warning_msg)
            return self

        seqinfos = self.seqInfo
        seqlengths = seqinfos.seqlengths
        isCircular = seqinfos.isCircular

        if seqlengths is None:
            print(warning_msg)
            return self

        if isCircular is None:
            print("considering all sequences as non-circular...")

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

        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        if ignoreStrand:
            obj["strand"] = ["*"] * len(self.column("seqnames"))

        df_gr = pd.DataFrame(obj)
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])

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
            withRevMap=withRevMap, minGapwidth=10000000, ignoreStrand=ignoreStrand
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
    ) -> "GenomicRanges":
        """Identify gaps in genomic positions for each distinct chromosome and strand combination.

        Args:
            start (int, optional): restrict chromosome start position. Defaults to 1.
            end (Optional[MutableMapping[str, int]], optional): restrict end position for each chromosome. Defaults to None.
                If None, this tried to use the SeqInfo object if available.

        Returns:
            GenomicRanges: A new GenomicRanges containing gap regions across chromosome and strand
        """

        # seqlengths = self._computeSeqLengths()
        seqlengths = self.seqlengths

        if end is None:
            end = seqlengths

        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        df_gr = pd.DataFrame(obj)
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])
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
    ) -> "GenomicRanges":
        """Calculate disjoint genomic positions for each distinct chromosome and strand combination.

        Args:
            withRevMap (bool, optional):return map of indices back to original object? . Defaults to False.
            ignoreStrand (bool, optional): ignore strand?. Defaults to False.

        Returns:
            GenomicRanges: A new GenomicRanges containing disjoint ranges across chromosome and strand
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

        df_gr = pd.DataFrame(obj)
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])
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
        a_obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        a_df_gr = pd.DataFrame(a_obj)

        b_obj = {
            "seqnames": other.column("seqnames"),
            "starts": other.column("starts"),
            "ends": other.column("ends"),
            "strand": other.column("strand"),
            "index": range(len(other.column("seqnames"))),
        }

        b_df_gr = pd.DataFrame(b_obj)

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

        if len(union_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = pd.DataFrame.from_records(union_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return GenomicRanges.fromPandas(final_df)

    def intersect(self, x: "GenomicRanges") -> "GenomicRanges":
        """Find intersection of genomic intervals with `other`

        Args:
            other (GenomicRanges): the other GenomicRanges object

        Returns:
            GenomicRanges: a new GenomicRanges object with 
                the intervals
        """
        a_obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        a_df_gr = pd.DataFrame(a_obj)

        b_obj = {
            "seqnames": x.column("seqnames"),
            "starts": x.column("starts"),
            "ends": x.column("ends"),
            "strand": x.column("strand"),
            "index": range(len(x.column("seqnames"))),
        }

        b_df_gr = pd.DataFrame(b_obj)

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

    def setdiff(self, x: "GenomicRanges") -> "GenomicRanges":
        """Find set difference of genomic intervals with `other`

        Args:
            other (GenomicRanges): the other GenomicRanges object

        Returns:
            GenomicRanges: a new GenomicRanges object with 
                the diff intervals
        """
        a_obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        a_df_gr = pd.DataFrame(a_obj)

        b_obj = {
            "seqnames": x.column("seqnames"),
            "starts": x.column("starts"),
            "ends": x.column("ends"),
            "strand": x.column("strand"),
            "index": range(len(x.column("seqnames"))),
        }

        b_df_gr = pd.DataFrame(b_obj)

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

    def binnedAverage(bins, numvar: str, varname: str):
        pass

    def subtract(
        self, x: "GenomicRanges", minoverlap: int = 1, ignoreStrand: bool = False
    ):
        pass

    # integer range methods
    def coverage(
        self, method: str, shift: int = 0, width: Optional[int] = None, weight: int = 1
    ):
        pass

    # search based methods
    def findOverlaps(
        self,
        query: "GenomicRanges",
        rtype: str,
        select: str,
        numResults: int = 1,
        mapGap: int = -1,
        minOverlap: int = 0,
        ignoreStrand: bool = False,
    ):
        pass

    def countOverlaps(
        self,
        query: "GenomicRanges",
        rtype: str,
        select: str,
        numResults: int = 1,
        mapGap: int = -1,
        minOverlap: int = 0,
        ignoreStrand: bool = False,
    ):
        pass

    def subsetByOverlaps(
        self,
        query: "GenomicRanges",
        rtype: str,
        select: str,
        numResults: int = 1,
        mapGap: int = -1,
        minOverlap: int = 0,
        ignoreStrand: bool = False,
    ):
        pass

    def nearest(
        self,
        query: "GenomicRanges",
        select: str = None,
        k: int = 1,
        ignoreStrand: bool = False,
    ) -> List[Optional[List[int]]]:
        """Find nearest positions that overlap with the each genomics interval in `x`. 
            For each interval in `x`, returns list of indices that match.

        Args:
            x (GenomicRanges): input GenomicRanges to find nearest positions
            k (int, optional): find k nearest positions. Defaults to 1.

        Returns:
            List[Optional[List[int]]]: List of possible hit indices for each interval in `x`
        """
        hits = []
        for _, row in query:
            grhits = self._index[row["seqnames"]].find_overlap(
                row["starts"], row["ends"]
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

    def precede(self, query: "GenomicRanges", select: str):
        pass

    def follow(self, query: "GenomicRanges", select: str):
        pass

    def distanceToNearest(self, query: "GenomicRanges"):
        pass

    # compare and order methods

    def duplicated(
        self,
        query: "GenomicRanges",
        method: str = "auto",
        incomparables: bool = False,
        fromLast: bool = False,
        nmax: Optional[int] = None,
    ):
        pass

    def match(query, nomatch: int = -1000000, incomparables: Optional[bool] = None):
        pass

    def isUnsorted(self, naRM=False, strictly=False, ignoreStrand=False):
        pass

    def order(naLast=True, decreasing=False, method: str = "auto"):
        pass

    def sort(self, by: str, decreasing: bool = False, ignoreStrand: bool = False):
        pass

    def rank(self, tiesMethod: str, naLast=False, method: str = "auto"):
        pass

    def pcompare(self, query: "GenomicRanges"):
        pass

    # windowing functions
    def tile(self, n: int, width: int):
        pass

    def slidingWindows(self, width: int, step: int = 1):
        pass

    def tileGenome(self, n: int, tilewidth: int, cutLastTileInChrom: bool = False):
        pass

    # summary methods
    def table(self, colName: Optional[str] = None):
        pass

    def sum(self, colName: str):
        pass

    def summary(self, colName: str):
        pass

    def sample(self, k: int = 5) -> "GenomicRanges":
        pass

    # misc
    def invertStrand(self) -> "GenomicRanges":
        pass

    def combine(self, **kwargs) -> "GenomicRanges":
        pass

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
