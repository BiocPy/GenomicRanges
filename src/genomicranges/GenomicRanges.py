import math
import random
from collections import OrderedDict
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Literal,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Union,
)
from warnings import warn

from biocframe import BiocFrame
from biocframe.types import SlicerArgTypes
from numpy import concatenate, count_nonzero, ndarray, sum, zeros
from pandas import DataFrame, concat, isna
from prettytable import PrettyTable

from .io import from_pandas
from .SeqInfo import SeqInfo
from .utils import (
    OVERLAP_QUERY_TYPES,
    calc_row_gapwidth,
    compute_mean,
    create_np_interval_vector,
    find_diff,
    find_disjoin,
    find_gaps,
    find_intersect,
    find_nearest,
    find_overlaps,
    find_union,
    slide_intervals,
    split_intervals,
)

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRanges(BiocFrame):
    """``GenomicRanges`` provides functionality to represent and operate over genomic regions and annotations.

    **Note: Intervals are inclusive on both ends and start at 1.**

    Additionally, ``GenomicRanges`` may also contain `Sequence Information` (checkout
    :py:class:`~genomicranges.SeqInfo.SeqInfo`) as part of its metadata. It contains for each
    sequence name (or chromosome) in the gene model, its length. Additionally, (checkout
    :py:class:`~genomicranges.SeqInfo.SeqInfo`) might also contain metadata about the
    genome, e.g. if its circular (`is_circular`) or not.

    Note: The documentation for some of the methods come from the
    `GenomicRanges R/Bioconductor package <https://github.com/Bioconductor/GenomicRanges>`_.

    Typical usage example:

    To construct a **GenomicRanges** object, simply pass in the column representation as a
    dictionary. This dictionary must contain "seqnames", "starts", "ends" columns and optionally
    specify "strand". If "strand" column is not provided, "*" is used as the default value for
    each genomic interval.

    .. code-block:: python

        gr = GenomicRanges(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [100, 115, 119],
                "ends": [103, 116, 120],
            }
        )

    Alternatively, you may also convert a :py:class:`~pandas.DataFrame` to ``GenomicRanges``.

    .. code-block:: python

        df = pd.DataFrame(
            {
                "seqnames": ["chr1", "chr2", "chr3"],
                "starts": [100, 115, 119],
                "ends": [103, 116, 120],
            }
        )

        gr = genomicranges.from_pandas(df)

    All columns other than "seqnames", "starts", "ends" and "strand" are considered
    metadata columns and can be accessed by
    :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.mcols`.

    .. code-block:: python

        gr.mcols()

    or slice the object

    .. code-block:: python

        sliced_gr = gr[1:2, [True, False, False]]

    Attributes:
        data (Mapping[str, Union[List[Any], Mapping]]):
            Columns as dictionary, Must contain "seqnames", "starts", "ends" and "strand" columns.
        number_of_rows (int, optional): Number of genomic intervals (or rows). Defaults to None.
        row_names (Sequence[str], optional): Row index. Defaults to None.
        column_names (Sequence[str], optional): Column names, automatically inferred from `data`.
            Defaults to None.
        metadata (Mapping, optional): Additional metadata. Defaults to None.
    """

    required_columns = ["seqnames", "starts", "ends", "strand"]

    def __init__(
        self,
        data: Mapping[str, Union[Sequence, Mapping]],
        number_of_rows: Optional[int] = None,
        row_names: Optional[Sequence] = None,
        column_names: Optional[Sequence] = None,
        metadata: Optional[Mapping] = None,
    ) -> None:
        """Initialize a `GenomicRanges` object."""
        super().__init__(data, number_of_rows, row_names, column_names, metadata)

    def _validate(self):
        """Internal function to validate ``GenomicRanges``."""
        if "strand" not in self._data:
            self._data["strand"] = ["*"] * len(self._data["starts"])

            if self.column_names is not None:
                self.column_names.append("strand")

        super()._validate()
        self._validate_ranges()

    def _validate_ranges(self):
        """Internal function to validate all columns of ``GenomicRanges``.

        Raises:
            ValueError: If missing required columns.
        """
        missing = list(set(self.required_columns).difference(set(self.column_names)))

        if len(missing) > 0:
            raise ValueError(
                f"`data` must contain {', '.join(self.required_columns)}."
                f"missing {missing} column{'s' if len(missing) > 1 else ''}"
            )

    @property
    def seqnames(self) -> List[str]:
        """Get sequence or chromosome names.

        Returns:
            List[str]: List of all chromosome names.
        """
        return self.column("seqnames")

    @property
    def start(self) -> List[int]:
        """Get sequence or chromosome start positions.

        Returns:
            List[int]: List of all chromosome start positions.
        """
        return self.column("starts")

    @property
    def end(self) -> List[int]:
        """Get sequence or chromosome end positions.

        Returns:
            List[int]: List of all chromosome end positions.
        """
        return self.column("ends")

    def ranges(
        self, ignore_strand: bool = False, return_type: Optional[Callable] = None
    ) -> Union[DataFrame, Dict, "GenomicRanges", Any]:
        """Get genomic positions.

        Args:
            ignore_strand (bool): Whether to ignore strands. Defaults to False.
            return_type (Callable, optional): Format to return genomic positions.
                Defaults to a dictionary representation, supports `DataFrame`
                or any callable representation that takes a dictionary as an input.

        Raises:
            ValueError: If ``return_type`` is not supported.

        Returns:
            Union[DataFrame, MutableMapping, "GenomicRanges", Any]: Genomic regions.
        """

        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
        }

        if ignore_strand:
            obj["strand"] = ["*"] * len(obj["seqnames"])

        if return_type is None:
            return obj
        else:
            try:
                return return_type(obj)
            except Exception as e:
                raise ValueError(f"{return_type} is not supported, {str(e)}")

    @property
    def strand(self) -> List[str]:
        """Get strand information.

        If ``strand`` is not provided, we use '*' as the default value for each interval.

        Returns:
            List[str]: Strand across all positions.
        """
        return self.column("strand")

    @property
    def width(self) -> List[int]:
        """Get widths of each interval.

        Returns:
            List[int]: Widths of each interval.
        """

        widths = []

        for _, row in self:
            widths.append(row["ends"] - row["starts"])

        return widths

    @property
    def seq_info(self) -> Optional[SeqInfo]:
        """Get sequence information, if available.

        Returns:
            (SeqInfo, optional): Sequence information.
        """

        if self.metadata and "seq_info" in self.metadata:
            return self.metadata["seq_info"]

        return None

    @seq_info.setter
    def seq_info(self, seq_info: Optional[SeqInfo]):
        """Set sequence information.

        Args:
            (SeqInfo, optional): Sequence information.

        Raises:
            ValueError: If `seq_info` is not a `SeqInfo` class.
        """

        if seq_info is not None:
            if not isinstance(seq_info, SeqInfo):
                raise ValueError("seq_info is not a `SeqInfo` class.")

            if self.metadata is None:
                self.metadata = {}

        self.metadata["seq_info"] = seq_info

    @property
    def seqlengths(self) -> Optional[Dict[str, int]]:
        """Get length of each chromosome, if available.

        Returns:
            (Dict[str, int], optional): A dictionary where keys are chromosome names and values
            specify lengths of each chromsome.
        """

        if self.metadata is not None and "seq_info" in self.metadata:
            return self.metadata["seq_info"].seqlengths

        return None

    @property
    def score(self) -> Optional[Sequence[Union[int, float]]]:
        """Get "score" column (if available) for each genomic interval.

        Returns:
            (Sequence[Union[int, float]], optional): Score column.
        """

        if "score" in self.column_names:
            return self.data["score"]

        return None

    @score.setter
    def score(self, score: Sequence[Union[int, float]]):
        """Set score for each position.

        Args:
            score (Sequence[Union[int, float]]): Score values to set.

        Raises:
            ValueError: If length of ``score`` does not
                match the number of intervals in the object.
            TypeError: If `score` is not a list.
        """

        if not isinstance(score, list):
            raise TypeError("`score` must be a list!")

        if len(score) != self.shape[0]:
            raise ValueError(
                "Length of `score` must be the same as number of intervals"
                f"must be {self.shape[0]}, but provided {len(score)}."
            )

        self["score"] = score

    @property
    def is_circular(self) -> Optional[Dict[str, bool]]:
        """Whether the sequences/chromosomes are circular (only if available).

        Returns:
            (Dict[str, bool], optional): A dictionary with chromosome names as keys and a
            boolean value indicating if its circular or not.
        """

        if self.metadata is not None and "seqInfo" in self.metadata:
            return self.metadata["seqInfo"].is_circular

        return None

    @property
    def genome(self) -> Optional[str]:
        """Get genome information, if available.

        Returns:
            (str, optional): Genome information if available else None.
        """

        if self._metadata and "seqInfo" in self._metadata:
            return self._metadata["seqInfo"].genome

        return None

    def granges(self) -> "GenomicRanges":
        """Creates a new ``GenomicRanges`` object with only ranges (`seqnames`, `starts, `ends` and `strand`).

        Returns:
            GenomicRanges: A new `GenomicRanges` with only ranges.
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
        self, return_type: Optional[Callable] = None
    ) -> Union[DataFrame, MutableMapping, Any]:
        """Get metadata across all genomic intervals.

        All columns other than "seqnames", "starts", "ends" and `"strand"
        are considered metadata for each interval.

        Args:
            return_type (Callable, optional): Format to return metadata.
                Defaults to dictionary representation, supports `DataFrame`
                or any callable representation that takes a dictionary as an input.

        Raises:
            ValueError: If ``return_type`` is not supported.

        Returns:
            Union[DataFrame, MutableMapping, Any]: Metadata columns without
            genomic positions.
        """

        new_data = OrderedDict()
        for k in self.column_names:
            if k not in self.required_columns:
                new_data[k] = self.column(k)

        if return_type is None:
            return new_data
        else:
            try:
                return return_type(new_data)
            except Exception as e:
                raise ValueError(f"{return_type} not supported, {str(e)}")

    def __repr__(self) -> str:
        table = PrettyTable(padding_width=2)
        table.field_names = [str(col) for col in self.column_names]

        _rows = []
        rows_to_show = 2
        _top = self.shape[0]
        if _top > rows_to_show:
            _top = rows_to_show

        # top three rows
        for r in range(_top):
            _row = self.row(r)
            vals = list(_row.values())
            res = [str(v) for v in vals]
            _rows.append(res)

        if self.shape[0] > 2 * rows_to_show:
            # add ...
            _rows.append(["..." for _ in range(len(self.column_names))])

        _last = self.shape[0] - rows_to_show
        if _last <= rows_to_show:
            _last = self.shape[0] - _top

        # last three rows
        for r in range(_last, len(self)):
            _row = self.row(r)
            vals = list(_row.values())
            res = [str(v) for v in vals]
            _rows.append(res)

        table.add_rows(_rows)

        pattern = (
            f"Class GenomicRanges with {self.dims[0]} intervals and "
            f"{self.dims[1] - 4} metadata columns \n"
            f"contains row names?: {self.row_names is not None} \n"
            f"{table.get_string()}"
        )

        return pattern

    # for documentation, otherwise serves no real use.
    def __getitem__(self, args: SlicerArgTypes) -> Union["GenomicRanges", dict, list]:
        """Subset the object.

        This operation returns a new object with the same type as caller.
        If you need to access specific rows or columns, use the
        :py:meth:`~biocframe.BiocFrame.BiocFrame.row` or
        :py:meth:`~biocframe.BiocFrame.BiocFrame.column`
        methods.

        Usage:

        .. code-block:: python

            # made up chromosome locations and ensembl ids.
            obj = {
                "ensembl": ["ENS00001", "ENS00002", "ENS00002"],
                "symbol": ["MAP1A", "BIN1", "ESR1"],
                "ranges": BiocFrame({
                    "chr": ["chr1", "chr2", "chr3"]
                    "start": [1000, 1100, 5000],
                    "end": [1100, 4000, 5500]
                ),
            }
            gr = GenomicRanges(obj)

            # different ways to slice.
            gr[0:2, 0:2]
            gr[[0,2], [True, False, False]]
            gr[<List of column names>]

        Args:
            args (SlicerArgTypes): A Tuple of slicer arguments to subset rows and
                columns. An element in ``args`` may be,

                - List of booleans, True to keep the row/column, False to remove.
                    The length of the boolean vector must be the same as number of rows/columns.

                - List of integer positions along rows/columns to keep.

                - A :py:class:`slice` object specifying the list of indices to keep.

                - A list of index names to keep. For rows, the object must contain unique
                    :py:attr:`~biocframe.BiocFrame.BiocFrame.row_names` and for columns must
                    contain unique :py:attr:`~biocframe.BiocFrame.BiocFrame.column_names`.

                - An integer to subset either a single row or column index.
                    Alternatively, you might want to use
                    :py:meth:`~biocframe.BiocFrame.BiocFrame.row` or
                    :py:meth:`~biocframe.BiocFrame.BiocFrame.column` methods.

                - A string to subset either a single row or column by label.
                    Alternatively, you might want to use
                    :py:meth:`~biocframe.BiocFrame.BiocFrame.row` or
                    :py:meth:`~biocframe.BiocFrame.BiocFrame.column` methods.

        Raises:
            ValueError: Too many slices provided.
            TypeError: If provided ``args`` are not an expected type.

        Returns:
            Union["GenomicRanges", dict, list]:
            - If a single row is sliced, returns a :py:class:`dict`.
            - If a single column is sliced, returns a :py:class:`list`.
            - For all other scenarios, returns the same type as caller with the subsetted rows and columns.
        """
        return super().__getitem__(args)

    # intra-range methods
    def flank(
        self,
        width: int,
        start: bool = True,
        both: bool = False,
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
        """Compute flanking ranges for each range. The logic for this comes from the `R/GenomicRanges` & `IRanges`
        packages.

        If ``start`` is ``True`` for a given range, the flanking occurs at the `start`,
        otherwise the `end`.
        The `widths` of the flanks are given by the ``width`` parameter.

        ``width`` can be negative, in which case the flanking region is
        reversed so that it represents a prefix or suffix of the range.

        Usage:

            `gr.flank(3, True)`, where "x" indicates a range in ``gr`` and "-" indicates the
            resulting flanking region:
                ---xxxxxxx
            If ``start`` were ``False``, the range in ``gr`` becomes
                xxxxxxx---
            For negative width, i.e. `gr.flank(x, -3, FALSE)`, where "*" indicates the overlap
            between "x" and the result:
                xxxx***
            If ``both`` is ``True``, then, for all ranges in "x", the flanking regions are
            extended into (or out of, if ``width`` is negative) the range, so that the result
            straddles the given endpoint and has twice the width given by width.

            This is illustrated below for `gr.flank(3, both=TRUE)`:
                ---***xxxx

        Args:
            width (int): Width to flank by. May be negative.
            start (bool, optional): Whether to only flank starts. Defaults to True.
            both (bool, optional): Whether to flank both starts and ends. Defaults to False.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            GenomicRanges: A new :py:class:`~.GenomicRanges` object with the flanked ranges.
        """
        new_starts = []
        new_ends = []

        all_starts = self.column("starts")
        all_ends = self.column("ends")
        all_strands = self.column("strand")

        # figure out which position to pin, start or end?
        start_flags = [start] * len(all_strands)
        if not ignore_strand:
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
                    tstart = all_starts[idx] if sf else all_ends[idx] + width + 1

            new_starts.append(tstart)
            new_ends.append(tstart + (width * (2 if both else 1) - 1))

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def resize(
        self,
        width: int,
        fix: Literal["start", "end", "center"] = "start",
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
        """Resize ranges to the specified ``width`` where either the ``start``, ``end``, or ``center`` is used as an
        anchor.

        Args:
            width (int): Width to resize, cannot be negative!
            fix (Literal["start", "end", "center"], optional): Fix positions by
                "start", "end", "center". Defaults to "start".
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Raises:
            ValueError: If parameter ``fix`` is neither `start`, `end` or `center`.
            ValueError: If ``width`` is negative.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with the resized ranges.
        """

        if width < 0:
            raise ValueError("`width` cannot be negative!")

        if fix not in ["start", "end", "center"]:
            raise ValueError(
                f"`fix` must be either 'start', 'end' or 'center', provided {fix}"
            )

        new_starts = []
        new_ends = []

        for idx, row in self:
            ts = None
            te = None

            if ignore_strand is True or row["strand"] != "-":
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
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def shift(self, shift: int = 0) -> "GenomicRanges":
        """Shift all intervals by parameter ``shift``.

        Args:
            shift (int, optional): Shift interval. Defaults to 0.
                If shift is 0, the current object is returned.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with the shifted ranges.
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
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def promoters(self, upstream: int = 2000, downstream: int = 200) -> "GenomicRanges":
        """Extend intervals to promoter regions.

        Generates promoter ranges relative to the transcription start site (TSS),
        where TSS is start(x). The promoter range is expanded around the TSS
        according to the upstream and downstream arguments. Upstream represents
        the number of nucleotides in the 5' direction and downstream the number
        in the 3' direction. The full range is defined as, (`start(x) - upstream`)
        to (`start(x) + downstream - 1`).

        Args:
            upstream (int, optional): Number of positions to extend in the 5' direction.
                Defaults to 2000.
            downstream (int, optional): Number of positions to extend in the 3' direction.
                Defaults to 200.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with ranges replaced
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
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def restrict(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        keep_all_ranges: bool = False,
    ) -> "GenomicRanges":
        """Restrict ranges to a given start and end positions.

        Args:
            start (int, optional): Start position. Defaults to None.
            end (int, optional): End position. Defaults to None.
            keep_all_ranges (bool, optional): Whether to keep intervals that do
                not overlap with start and end. Defaults to False.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with restricted ranges.
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

        new_row_names = self.row_names.copy()

        if not keep_all_ranges:
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

            for col in self.column_names:
                new_data[col] = [new_data[col][x] for x in keepers]

            if self.row_names:
                new_row_names = [new_row_names[x] for x in keepers]

        return GenomicRanges(
            new_data,
            row_names=new_row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def trim(self) -> "GenomicRanges":
        """Trim sequences outside of bounds for non-circular chromosomes.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with trimmed ranges.
        """

        if self.seq_info is None:
            raise ValueError("Cannot trim ranges. `seqinfo` is not available.")

        if self.metadata is None:
            raise ValueError("Cannot trim ranges. `seqinfo` is not available.")

        if "seq_info" in self.metadata and self.metadata["seq_info"] is None:
            raise ValueError("Cannot trim ranges. `seqinfo` is not available.")

        seqinfos = self.seq_info
        seqlengths = seqinfos.seqlengths
        is_circular = seqinfos.is_circular

        if seqlengths is None:
            raise ValueError("Cannot trim ranges. `seqlengths` is not available.")

        if is_circular is None:
            warn("considering all sequences as non-circular...")

        all_chrs = self.column("seqnames")
        all_ends = self.column("ends")

        new_data = self._data.copy()
        new_row_names = self.row_names.copy()

        keepers = []
        for idx in range(len(all_chrs)):
            keep = True
            t_chr = all_chrs[idx]
            if (
                is_circular is not None
                and is_circular[t_chr] is False
                and all_ends[idx] > seqlengths[t_chr]
            ):
                keep = False

            if keep:
                keepers.append(idx)

        for col in self.column_names:
            new_data[col] = [new_data[col][x] for x in keepers]

        if self.row_names:
            new_row_names = [new_row_names[x] for x in keepers]

        return GenomicRanges(
            new_data,
            row_names=new_row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    # TODO: needs checks when relative - {start, width and end} do not agree
    def narrow(
        self,
        start: Optional[int] = None,
        width: Optional[int] = None,
        end: Optional[int] = None,
    ) -> "GenomicRanges":
        """Narrow genomic positions by provided ``start``, ``width`` and ``end`` parameters.

        Important: these parameters are relative shift in positions for each range.

        Args:
            start (int, optional): Relative start position. Defaults to None.
            width (int, optional): Relative end position. Defaults to None.
            end (int, optional): Relative width of the interval. Defaults to None.

        Raises:
            ValueError: If `width` is provided, either `start` or `end` must be provided.
            ValueError: Provide two of the three parameters - `start`, `end` and `width`
                but not all.

        Returns:
            GenomicRanges:  A new `GenomicRanges` object with narrow ranges.
        """
        if start is not None and end is not None and width is not None:
            raise ValueError(
                "Only provide two of the three parameters - `start`, "
                "`end` and `width` but not all!"
            )

        if width is not None:
            if start is None and end is None:
                raise ValueError(
                    "If width is provided, either start or end must be provided."
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
                    "If `width` is provided, either `start` or `end` must be provided."
                )

        new_data = self._data.copy()
        new_data["starts"] = new_starts
        new_data["ends"] = new_ends

        return GenomicRanges(
            new_data,
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def _calc_gap_widths(self, ignore_strand: bool = False) -> List[int]:
        """Internal method to calculate gap widths.

        Args:
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            List[int]: Gap widths for each range.
        """
        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": [i for i in range(len(self.column("seqnames")))],
        }

        df_gr = DataFrame(obj)
        df_gr = df_gr.sort_values(by=["seqnames", "strand", "starts", "ends"])

        if ignore_strand:
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
    # TODO: this is a very ineffecient implementation, can do a better job later.
    def reduce(
        self,
        with_reverse_map: bool = False,
        min_gap_width: int = 1,
        ignore_strand: bool = False,
    ) -> "GenomicRanges":
        """Reduce orders the ranges, then merges overlapping or adjacent ranges.

        Args:
            with_reverse_map (bool, optional): Whether to return map of indices back to
                original object. Defaults to False.
            min_gap_width (int, optional): Ranges separated by a gap of
                at least ``min_gap_width`` positions are not merged. Defaults to 1.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with reduced intervals.
        """

        df_gr = self._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)

        df_gr["gapwidths"] = self._calc_gap_widths(ignore_strand=ignore_strand)
        df_gr["gapwidth_flag"] = [
            (True if (isna(x) or (x == 0)) else x < min_gap_width)
            for x in df_gr["gapwidths"]
        ]

        gaps_to_merge = df_gr[df_gr["gapwidth_flag"] == True]  # noqa: E712

        gaps_merged = gaps_to_merge.groupby(
            ["seqnames", "strand", "gapwidth_flag"], sort=False
        ).agg(starts=("starts", min), ends=("ends", max), revmap=("index", list))

        gaps_merged = gaps_merged.reset_index()

        gaps_not_merged = df_gr[df_gr["gapwidth_flag"] == False]  # noqa: E712
        gaps_not_merged["revmap"] = gaps_not_merged["index"].apply(lambda x: [x])
        gaps_not_merged = gaps_not_merged[gaps_merged.columns]

        finale = concat([gaps_merged, gaps_not_merged])

        columns_to_keep = ["seqnames", "strand", "starts", "ends"]

        if with_reverse_map:
            columns_to_keep.append("revmap")

        finale = finale[columns_to_keep].sort_values(
            ["seqnames", "strand", "starts", "ends"]
        )

        return from_pandas(finale)

    def range(
        self, with_reverse_map: bool = False, ignore_strand: bool = False
    ) -> "GenomicRanges":
        """Calculate ranges for each chromosome. (minimum of all starts, maximum of all ends) in the object.

        Technically its same as :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.reduce`
        with a ridiculously high ``min_gap_width``.

        Args:
            with_reverse_map (bool, optional): return map of indices back to
                original object?. Defaults to False.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            GenomicRanges: a new `GenomicRanges` object with the new ranges.
        """
        return self.reduce(
            with_reverse_map=with_reverse_map,
            min_gap_width=100000000000,
            ignore_strand=ignore_strand,
        )

    def compute_seqlengths(self) -> Dict[str, int]:
        """Get seqlengths either from the :py:class:`~genomicranges.SeqInfo.SeqInfo` object or computes one from the
        current ranges.

        Note: If computed, they are specific to this `GenomicRanges` object
        and may not represent the seqlenths of the genome.

        Returns:
            Dict[str, int]: A dict with chromosome names as keys and their lengths as values.
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
        """Identify gaps in genomic positions for each distinct `seqname` (chromosome) and `strand` combination.

        Args:
            start (int, optional): Restrict chromosome start position. Defaults to 1.
            end (MutableMapping[str, int], optional): Restrict end
                position for each chromosome. Defaults to None. If None, it uses the
                :py:class:`~genomicranges.SeqInfo.SeqInfo` object if available.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` containing
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
            groups = df_gr.groupby("seqnames")
            for name, group in groups:
                strands = group["strand"].unique()
                missing_strands = list(set(["*", "+", "-"]).difference(set(strands)))
                for ms in missing_strands:
                    gap_intervals.append((name, ms, start, end[name]))

        if len(gap_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        final_df = DataFrame.from_records(gap_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def disjoin(
        self, with_reverse_map: bool = False, ignore_strand: bool = False
    ) -> Optional["GenomicRanges"]:
        """Calculate disjoint genomic positions for each distinct `seqname` (chromosome) and `strand` combination.

        Args:
            with_reverse_map (bool, optional): Whether to return map of indices back to
                original object. Defaults to False.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` containing
            disjoint ranges across chromosome and strand.
        """

        df_gr = self._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)
        groups = df_gr.groupby(["seqnames", "strand"])

        disjoin_intervals = []
        for name, group in groups:
            group["iindex"] = range(len(group))

            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            tdisjoin = find_disjoin(all_intvals, with_reverse_map=with_reverse_map)

            for td in tdisjoin:
                td_res = (name[0], name[1], td[0], td[1])
                if with_reverse_map:
                    td_res.append(td[2])
                disjoin_intervals.append(td_res)

        if len(disjoin_intervals) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends"]
        if with_reverse_map:
            columns.append("revmap")

        final_df = DataFrame.from_records(disjoin_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def is_disjoint(self, ignore_strand: bool = False) -> bool:
        """Whether all ranges (for each chromosome, strand) are disjoint.

        Args:
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            bool: True if the intervals are disjoint.
        """
        new_gr = self.disjoin(ignore_strand=ignore_strand)

        return new_gr.shape != self.shape

    def disjoint_bins(self, ignore_strand: bool = False):
        raise NotImplementedError

    # set operations
    def union(self, other: "GenomicRanges") -> "GenomicRanges":
        """Find union of genomic intervals with `other`.

        Args:
            other (GenomicRanges): The other `GenomicRanges` object.

        Raises:
            TypeError: If ``other`` is not a `GenomicRanges`.

        Returns:
            GenomicRanges: A new `GenomicRanges` object with all ranges.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("`other` is not a `GenomicRanges` object")

        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        df_gr = concat([a_df_gr, b_df_gr])
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
        final_df = DataFrame.from_records(union_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def intersect(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find intersection of genomic intervals with `other`.

        Args:
            other (GenomicRanges): The other `GenomicRanges` object.

        Raises:
            TypeError: If ``other`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` object with intersecting ranges.
            If there are no overlapping intervals, returns None.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("`other` is not a `GenomicRanges` object")

        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        df_gr = concat([a_df_gr, b_df_gr])
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
        final_df = DataFrame.from_records(intersect_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def setdiff(self, other: "GenomicRanges") -> Optional["GenomicRanges"]:
        """Find set difference of genomic intervals with `other`.

        Args:
            other (GenomicRanges): The other `GenomicRanges` object.

        Raises:
            TypeError: If ``other`` is not of type `GenomicRanges`.


        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` object with the diff ranges.
            If there are no diff intervals, returns None.
        """

        if not isinstance(other, GenomicRanges):
            raise TypeError("`other` is not a `GenomicRanges` object")

        a_df_gr = self._generic_pandas_ranges(sort=False)
        b_df_gr = other._generic_pandas_ranges(sort=False)

        only_seqnames = concat(
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
        final_df = DataFrame.from_records(diff_ints, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def binned_average(
        self, scorename: str, bins: "GenomicRanges", outname: str
    ) -> "GenomicRanges":
        """Calculate average for a column across all bins, then set a column called ``outname`` with those values.

        Args:
            scorename (str): Score column to compute averages on.
            bins (GenomicRanges): Bins you want to use.
            outname (str): New column name to add to the object.

        Raises:
            ValueError: If ``scorename`` column does not exist.
            Exception: ``scorename`` is not all ints or floats.
            TypeError: If ``bins`` is not of type `GenomicRanges`.

        Returns:
            GenomicRanges: A new `GenomicRanges` similar to bins,
            with a new column containing the averages.
        """

        if not isinstance(bins, GenomicRanges):
            raise TypeError("`bins` is not a `GenomicRanges` object.")

        if scorename not in self.column_names:
            raise ValueError(f"'{scorename}' is not a valid column name")

        values = self.column(scorename)

        if not all((isinstance(x, int) or isinstance(x, float)) for x in values):
            raise Exception(f"'{scorename}' values must be either ints or floats.")

        df_gr = self._generic_pandas_ranges(sort=True)
        df_gr["values"] = values

        tgt_gr = bins._generic_pandas_ranges(ignore_strand=True, sort=True)
        tgt_groups = tgt_gr.groupby("seqnames")

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
                vec_mean = sum(vec) / count_nonzero(vec)

                result.append((name, tint["starts"], tint["ends"], "*", vec_mean))

        columns = ["seqnames", "starts", "ends", "strand", outname]
        final_df = DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def subtract(
        self, x: "GenomicRanges", min_overlap: int = 1, ignore_strand: bool = False
    ):
        raise NotImplementedError

    # integer range methods
    def coverage(
        self, shift: int = 0, width: Optional[int] = None, weight: int = 1
    ) -> Dict[str, ndarray]:
        """Calculate coverage for each chromosome, For each position, counts the number of ranges that cover it.

        Args:
            shift (int, optional): Shift all genomic positions. Defaults to 0.
            width (int, optional): Restrict the width of all
                chromosomes. Defaults to None.
            weight (int, optional): Weight to use. Defaults to 1.

        Returns:
            Dict[str, ndarray]:  A dictionary with chromosome names as keys and the
            coverage vector as value.
        """
        df_gr = self._generic_pandas_ranges(sort=True)
        groups = df_gr.groupby("seqnames")

        shift_arr = None
        if shift > 0:
            shift_arr = zeros(shift)

        result = {}
        for name, group in groups:
            all_intvals = [
                (x[0], x[1])
                for x in zip(group["starts"].to_list(), group["ends"].to_list())
            ]

            cov, _ = create_np_interval_vector(
                intervals=all_intvals, with_reverse_map=False
            )

            if shift > 0:
                cov = concatenate((shift_arr, cov))

            if weight > 0:
                cov = cov * weight

            if width is not None:
                cov = cov[:width]

            result[name] = cov

        return result

    # search based methods
    def find_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Find overlaps between ``subject`` (self) and a ``query`` `GenomicRanges` object.

        Args:
            query (GenomicRanges): Query `GenomicRanges`.
            query_type (Literal["any", "start", "end", "within"], optional): Overlap query type,
                must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".
            max_gap (int, optional): Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).
            min_overlap (int, optional): Minimum overlap with query. Defaults to 1.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` object
            with the same length as query, containing hits to overlapping indices.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("`query` is not a `GenomicRanges` object.")

        if query_type not in OVERLAP_QUERY_TYPES:
            raise ValueError(
                f"'{query_type}' must be one of {', '.join(OVERLAP_QUERY_TYPES)}."
            )

        df_gr = self._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)
        tgt_gr = query._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)

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
                max_gap=max_gap,
                min_overlap=min_overlap,
                query_type=query_type,
            )

            for th in thits:
                tindices = [src_intvals_map[i - 1] for i in th[2]]
                result.append((chrom, ustrand, th[0][0], th[0][1], tindices))

        if len(result) == 0:
            return None

        columns = ["seqnames", "strand", "starts", "ends", "hits"]
        final_df = DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def count_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> List[int]:
        """Count overlaps between ``subject`` (self) and a ``query`` `GenomicRanges` object.

        Args:
            query (GenomicRanges): Query `GenomicRanges`.
            query_type (Literal["any", "start", "end", "within"], optional): Overlap query type,
                must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".
            max_gap (int, optional): Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).
            min_overlap (int, optional): Minimum overlap with query. Defaults to 1.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            List[int]: Number of overlaps for each range in query.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("`query` is not a `GenomicRanges` object.")

        result = self.find_overlaps(
            query=query,
            query_type=query_type,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
        )

        if result is None:
            return None

        hits = result.column("hits")
        return [len(ht) for ht in hits]

    def subset_by_overlaps(
        self,
        query: "GenomicRanges",
        query_type: Literal["any", "start", "end", "within"] = "any",
        max_gap: int = -1,
        min_overlap: int = 1,
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Subset ``subject`` (self) with overlaps in ``query`` `GenomicRanges` object.

        Args:
            query (GenomicRanges): Query `GenomicRanges`.
            query_type (Literal["any", "start", "end", "within"], optional): Overlap query type,
                must be one of

                - "any": Any overlap is good
                - "start": Overlap at the beginning of the intervals
                - "end": Must overlap at the end of the intervals
                - "within": Fully contain the query interval

                Defaults to "any".
            max_gap (int, optional): Maximum gap allowed in the overlap.
                Defaults to -1 (no gap allowed).
            min_overlap (int, optional): Minimum overlap with query. Defaults to 1.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` object
            containing only subsets.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("`query` is not a `GenomicRanges` object.")

        result = self.find_overlaps(
            query=query,
            query_type=query_type,
            max_gap=max_gap,
            min_overlap=min_overlap,
            ignore_strand=ignore_strand,
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
        ignore_strand: bool = False,
        stepstart: int = 3,
        stepend: int = 3,
    ) -> Optional["GenomicRanges"]:
        """Internal function to search ``self`` and a ``query`` `GenomicRanges` object.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` object that has
            the same length as query but contains `hits` to
            `indices` and `distance`.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("`query` is not a `GenomicRanges` object.")

        subject_gr = self._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)
        query_gr = query._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)

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
        final_df = DataFrame.from_records(result, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def nearest(
        self,
        query: "GenomicRanges",
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions both upstream and downstream that overlap with each range in ``query``.

        Adds a new column ("hits") to ``query`` with the nearest matches.

        Args:
            query (GenomicRanges): Query `GenomicRanges` to find nearest positions.
            ignore_strand (bool, optional): Whether to ignore strand. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` containing list of
            possible "hits" indices for each range in ``query``.
        """
        return self._generic_search(query=query, ignore_strand=ignore_strand)

    def precede(
        self,
        query: "GenomicRanges",
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with each range in ``query``.

        Adds a new column ("hits") to ``query`` with the nearest matches.

        Args:
            query (GenomicRanges): Query `GenomicRanges` to find nearest positions.
            ignore_strand (bool, optional): Whether to ignore strand. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` containing list of
            possible "hits" indices for each range in ``query``.
        """
        return self._generic_search(query=query, ignore_strand=ignore_strand, stepend=0)

    def follow(
        self,
        query: "GenomicRanges",
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only upstream that overlap with each range in ``query``.

        Adds a new column ("hits") to ``query`` with the nearest matches.

        Args:
            query (GenomicRanges): Query `GenomicRanges` to find nearest positions.
            ignore_strand (bool, optional): Whether to ignore strand. Defaults to False.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): a new `GenomicRanges` containing list of
            possible "hits" indices for each range in ``query``.
        """
        return self._generic_search(
            query=query, ignore_strand=ignore_strand, stepstart=0
        )

    def distance_to_nearest(
        self,
        query: "GenomicRanges",
        ignore_strand: bool = False,
    ) -> Optional["GenomicRanges"]:
        """Search nearest positions only downstream that overlap with each range in ``query``.

        Adds a new column ("hits") to ``query`` with the nearest matches.

        Note: Technically same as :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.nearest`
        since we also return `distances`.

        Args:
            query (GenomicRanges): Query `GenomicRanges` to find nearest positions.
            ignore_strand (bool, optional): Whether to ignore strand. Defaults to False.

        Raises:
            TypeError: If query is not of type `GenomicRanges`.

        Returns:
            ("GenomicRanges", optional): A new `GenomicRanges` containing list of
            possible "hits" indices for each range in ``query``.
        """
        return self._generic_search(query=query, ignore_strand=ignore_strand)

    # compare and order methods
    def duplicated(
        self,
    ) -> List[bool]:
        """Element wise comparison to find duplicate ranges.

        Returns:
            List[bool]: True if duplicated else False.
        """
        df = self._generic_pandas_ranges(sort=False)
        return df.duplicated().to_list()

    def match(self, query: "GenomicRanges") -> List[Optional[int]]:
        """Element wise comparison to find exact match ranges.

        Args:
            query (GenomicRanges): Input `GenomicRanges` to search matches.

        Raises:
            TypeError: If ``query`` is not of type `GenomicRanges`.

        Returns:
            (List[int], optional): List contianing index positions if able to match else None.
        """

        if not isinstance(query, GenomicRanges):
            raise TypeError("`query` is not a `GenomicRanges` object.")

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

    def _generic_pandas_ranges(self, ignore_strand=False, sort=False) -> DataFrame:
        """Internal function to create a :py:class:`~pandas.DataFrame` from ranges.

        Args:
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.
            sort (bool, optional): Whether to sort by region. Defaults to False.

        Returns:
            DataFrame: a pandas `DataFrame` of the genomic ranges.
        """
        obj = {
            "seqnames": self.column("seqnames"),
            "starts": self.column("starts"),
            "ends": self.column("ends"),
            "strand": self.column("strand"),
            "index": range(len(self.column("seqnames"))),
        }

        if ignore_strand:
            obj["strand"] = ["*"] * len(self.column("seqnames"))

        df = DataFrame(obj)

        if sort:
            df = df.sort_values(["seqnames", "strand", "starts", "ends"])

        return df

    def _generic_order(self, ignore_strand=False) -> List[int]:
        """Internal function to compute order.

        Args:
            query (GenomicRanges): Query `GenomicRanges` to search matches.

        Returns:
            List[int, optional]: Sorted index order.
        """
        sorted = self._generic_pandas_ranges(ignore_strand=ignore_strand, sort=True)
        new_order = sorted["index"]

        return new_order

    def is_unsorted(self, ignore_strand=False) -> bool:
        """Whether the genomic positions are unsorted.

        Args:
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            bool: True if unsorted else False.
        """
        order = self._generic_order(ignore_strand=ignore_strand)
        diff = order.diff()

        if any(diff != 1):
            return True

        return False

    def order(self, decreasing=False) -> List[int]:
        """Get the order of indices for sorting.

        Args:
            decreasing (bool, optional): Whether to sort in descending order. Defaults to False.

        Returns:
            List[int]: List with order of indices.
        """
        order = self._generic_order()

        if decreasing:
            order = order[::-1]
        return order.to_list()

    def sort(
        self, decreasing: bool = False, ignore_strand: bool = False
    ) -> "GenomicRanges":
        """Sort the `GenomicRanges` object.

        Args:
            decreasing (bool, optional): Whether to sort in descending order. Defaults to False.
            ignore_strand (bool, optional): Whether to ignore strands. Defaults to False.

        Returns:
            "GenomicRanges": A new sorted `GenomicRanges` object.
        """
        order = self._generic_order(ignore_strand=ignore_strand)

        if decreasing:
            order = order[::-1]

        new_order = order.to_list()
        return self[new_order, :]

    def rank(self) -> List[int]:
        """Get rank of the `GenomicRanges` object.

        For each range identifies its position is a sorted order.

        Returns:
            List[int]: List of indices identifying rank.
        """
        order = self._generic_order().to_list()
        rank = [order.index(x) for x in range(len(order))]
        return rank

    # windowing functions
    def tile_by_range(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each range into chunks by ``n`` (number of intervals) or ``width`` (intervals with equal width).

        Note: Either ``n`` or ``width`` must be provided, but not both.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a gneomic into chunks.

        Args:
            n (int, optional): Number of intervals to split into.
                Defaults to None.
            width (int, optional): Width of each interval. Defaults to None.

        Raises:
            ValueError: If both ``n`` or ``width`` are provided.

        Returns:
            "GenomicRanges": A new `GenomicRanges` with the split ranges.
        """
        if n is not None and width is not None:
            raise ValueError("Either `n` or `width` must be provided but not both.")

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
        final_df = DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def tile(
        self, n: Optional[int] = None, width: Optional[int] = None
    ) -> "GenomicRanges":
        """Split each interval by ``n`` (number of sub intervals) or ``width`` (intervals with equal width).

        Note: Either ``n`` or ``width`` must be provided but not both.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a gneomic into chunks, or
        :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.tile_by_range`.

        Args:
            n (int, optional): Number of intervals to split into. Defaults to None.
            width (int, optional): Width of each interval. Defaults to None.

        Raises:
            ValueError: If both ``n`` and ``width`` are provided.

        Returns:
            "GenomicRanges": A new `GenomicRanges` with the split ranges.
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
        final_df = DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    def sliding_windows(self, width: int, step: int = 1) -> "GenomicRanges":
        """Slide along each range by ``width`` (intervals with equal ``width``) and ``step``.

        Also, checkout :py:func:`~genomicranges.io.tiling.tile_genome` for splitting
        a gneomic into chunks, or
        :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.tile_by_range`.

        Args:
            width (int, optional): Width of each interval. Defaults to None.
            step (int, optional): Step interval, Defaults to 1.

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
        final_df = DataFrame.from_records(all_intervals, columns=columns)
        final_df = final_df.sort_values(["seqnames", "strand", "starts", "ends"])
        return from_pandas(final_df)

    # misc methods
    def sample(self, k: int = 5) -> "GenomicRanges":
        """Randomly sample ``k`` intervals.

        Args:
            k (int, optional): Number of intervals to sample. Defaults to 5.

        Returns:
            GenomicRanges: A new `GenomicRanges` with randomly sampled ranges.
        """
        sample = random.sample(range(len(self)), k=k)
        return self[sample, :]

    def invert_strand(self) -> "GenomicRanges":
        """Invert strand information for each interval.

        Conversion map:
            - "+" map to "-"
            - "-" becomes "+"
            - "*" stays the same

        Returns:
            GenomicRanges: A new `GenomicRanges` object.
        """
        convertor = {"+": "-", "-": "+", "*": "*"}
        inverts = [convertor[idx] for idx in self.column("strand")]

        new_data = self._data.copy()
        new_data["strand"] = inverts

        return GenomicRanges(
            new_data,
            number_of_rows=self.shape[0],
            row_names=self.row_names,
            column_names=self.column_names,
            metadata=self.metadata,
        )

    def concat(self, *granges: "GenomicRanges") -> "GenomicRanges":
        """Row-wise concatenate multiple `GenomicRanges` objects.

        Args:
            granges (GenomicRanges): Objects to concatenate.

        Raises:
            TypeError: If any ``granges`` are not "GenomicRanges".

        Returns:
            GenomicRanges: A new concatenated `GenomicRanges` object.
        """
        all_granges = [isinstance(gr, GenomicRanges) for gr in granges]

        if not all(all_granges):
            raise TypeError("all provided objects are not GenomicRanges objects")

        all_columns = [gr.column_names for gr in granges]
        all_columns.append(self.column_names)
        all_unique_columns = list(
            set([item for sublist in all_columns for item in sublist])
        )

        new_data = OrderedDict()
        for col in all_unique_columns:
            if col not in new_data:
                new_data[col] = []

            for gr in granges:
                if col in gr.column_names:
                    new_data[col].extend(gr.column(col))
                else:
                    new_data[col].extend([None] * len(gr))

            if col in self.column_names:
                new_data[col].extend(self.column(col))
            else:
                new_data[col].extend([None] * len(self))

        return GenomicRanges(new_data, column_names=all_unique_columns)
