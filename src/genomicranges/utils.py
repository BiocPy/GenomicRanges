from random import random
from typing import Any, List, MutableMapping, Optional, Tuple, Union

from numpy import diff, ndarray, split, where, zeros

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def calc_row_gapwidth(
    a: MutableMapping[str, Any], b: MutableMapping[str, Any]
) -> Optional[int]:
    """Calculates gap width for two genomic positions from :py:class:`~genomicranges.GenomicRanges.GenomicRanges`.

    a,b must contain keys `seqnames`, `strand`, `starts` and `ends`.

    Args:
        a, b (MutableMapping[str, Any]): Genomic row with positions.

        Usually a result of :py:meth:`~genomicranges.GenomicRanges.GenomicRanges.row`.

    Returns:
        int: gapwidth,

        - ``nan`` if the sequences are not comparable
        - 0 if they overlap
        - a number if there is a gap.

        returns None if gap does not exist.
    """

    if a["seqnames"] != b["seqnames"] or a["strand"] != b["strand"] != "-":
        return float("nan")

    if (a["starts"] <= b["starts"] <= a["ends"]) or (
        a["starts"] <= b["ends"] <= a["ends"]
    ):
        return 0

    return b["starts"] - a["ends"]


def np_split_groups_ne(
    ary: ndarray, step: int = 0
) -> Tuple[List[ndarray], List[ndarray]]:
    """Split a :py:class:`~numpy.ndarray` by consecutive values.

    Args:
        ary (ndarray): A numpy array.
        step (int, optional): Value to split consecutive intervals by.

    Returns:
        Tuple(List[ndarray], List[ndarray]): Tuple with indices
        and their coverage.
    """
    return [
        (x[0], x[1])
        for x in zip(
            split(range(len(ary)), where(diff(ary) != step)[0] + 1),
            split(ary, where(diff(ary) != step)[0] + 1),
        )
    ]


def create_np_interval_vector(
    intervals: List[Tuple[int, int]],
    with_reverse_map: bool = False,
    force_size: Optional[int] = None,
    dont_sum: bool = False,
    value: int = 1,
) -> Tuple[ndarray, Optional[List]]:
    """Represent intervals/calculate coverage.

    Args:
        intervals (List[Tuple[int, int]]): Input interval vector.
        with_reverse_map (bool, optional): Return map of indices? Defaults to False.
        force_size (int, optional): Force size of the array.
        dont_sum (bool, optional): Do not sum. Defaults to False.
        value (int, optional): Default value to increment. Defaults to 1.

    Returns:
        Tuple[ndarray, Optional[List]]: A numpy array representing
        coverage from the intervals and optionally the index map.
    """
    if len(intervals) < 1:
        return intervals

    max_end = force_size
    if max_end is None:
        max_end = max([x[1] for x in intervals])
    cov = zeros(max_end)

    revmap = None
    if with_reverse_map:
        revmap = [[] for _ in range(max_end)]

    for idx in range(len(intervals)):
        i = intervals[idx]

        if dont_sum:
            cov[i[0] - 1 : i[1]] = value
        else:
            cov[i[0] - 1 : i[1]] += value

        if with_reverse_map:
            _ = [revmap[x].append(idx + 1) for x in range(i[0] - 1, i[1])]

    return cov, revmap


def find_unary_union(
    intervals: List[Tuple[int, int]], with_reverse_map: bool = False
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Union of all overlapping intervals.

    Args:
        intervals (List[Tuple[int, int]]): Input interval vector.
        with_reverse_map (bool, optional): Return map of indices. Defaults to False.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List
        containing tuples with (start, end, indices) for each contiguous region.
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, with_reverse_map=with_reverse_map
    )

    # splitting intervals by where gaps=0 gives us the union of all intervals
    groups = [
        (x[0], x[1])
        for x in zip(
            split(range(len(np_intvals)), where(np_intvals == 0)[0] + 1),
            split(np_intvals, where(np_intvals == 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if cov[0] != 0:
            start = min(idx) + 1
            end = max(idx) + 1 if (max(idx) == len(np_intvals) - 1) else max(idx)
            result = (
                start,
                end,
            )

            if with_reverse_map:
                result.append(
                    list(
                        set(
                            [
                                item
                                for sublist in revmap[min(idx) : max(idx) + 1]
                                for item in sublist
                            ]
                        )
                    )
                )

            final.append(result)

    return final


def find_gaps(
    intervals: List[Tuple[int, int]], start_limit: int = 1, end_limit: int = None
) -> List[Tuple[int, int]]:
    """Find gaps across all intervals.

    Args:
        intervals (List[Tuple[int, int]]): Input interval vector.
        start_limit (int, optional): Start for identifying gaps. Defaults to 1.
        end_limit (int, optional): End interval. Defaults to None.

    Returns:
        List[Tuple[int, int]]: List of tuples with
        (start, end) for each gap.
    """

    min_start = start_limit
    if min_start is None:
        min_start = min([x[0] for x in intervals])

    max_end = end_limit
    if max_end is None:
        max_end = max(x[1] for x in intervals)

    np_intvals, _ = create_np_interval_vector(
        intervals=intervals, with_reverse_map=False, force_size=max_end
    )

    # splitting intervals by diff gives us the gaps
    groups = [
        (x[0], x[1])
        for x in zip(
            split(range(len(np_intvals)), where(diff(np_intvals) != 0)[0] + 1),
            split(np_intvals, where(diff(np_intvals) != 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if cov[0] == 0:
            start = min(idx) + 1
            end = max(idx) + 1

            if start < min_start:
                start = min_start

            if end - start < 0 or start > end:
                continue

            result = (start, end)

            final.append(result)

    return final


def find_disjoin(
    intervals: List[Tuple[int, int]], with_reverse_map: bool = False
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Find disjoint sets across all intervals.

    Args:
        intervals (List[Tuple[int, int]]): Input interval vector.
        with_reverse_map (bool, optional): Return map of indices. Defaults to False.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List
        of tuples with (start, end).
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, with_reverse_map=with_reverse_map
    )

    # its isn't enough, so lets add some randomless to distinguish each interval
    for itval in intervals:
        np_intvals[itval[0] - 1 : itval[1]] += int(10 * random())

    min_start = min([x[0] for x in intervals])
    max_end = max(x[1] for x in intervals)

    # splitting intervals by diff gives us the gaps
    groups = [
        (x[0], x[1])
        for x in zip(
            split(range(len(np_intvals)), where(diff(np_intvals) != 0)[0] + 1),
            split(np_intvals, where(diff(np_intvals) != 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if cov[0] != 0:
            start = min(idx) + 1
            end = max(idx) + 1

            if start >= min_start and end <= max_end:
                result = (
                    start,
                    end,
                )

                if with_reverse_map:
                    result.append(
                        list(
                            set(
                                [
                                    item
                                    for sublist in revmap[min(idx) : max(idx) + 1]
                                    for item in sublist
                                ]
                            )
                        )
                    )

                final.append(result)

    return final


def find_unary_intersect(
    intervals: List[Tuple[int, int]], with_reverse_map: bool = False, threshold: int = 1
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Intersection of all intervals.

    Args:
        intervals (List[Tuple[int, int]]): Input interval vector.
        with_reverse_map (bool, optional): Return map of indices?. Defaults to False.
        threshold (int, optional): Threshold for cutoff. Defaults to 1.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List
        containing tuples with (start, end, indices) for each contiguous region.
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, with_reverse_map=with_reverse_map
    )

    min_start = min([x[0] for x in intervals])
    max_end = max(x[1] for x in intervals)

    # splitting intervals by where gaps=0 gives us the union of all intervals
    groups = [
        (x[0], x[1])
        for x in zip(
            split(range(len(np_intvals)), where(np_intvals == 0)[0] + 1),
            split(np_intvals, where(np_intvals == 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if cov[0] > threshold:
            start = min(idx) + 1
            end = max(idx) + 1 if (max(idx) == len(np_intvals) - 1) else max(idx)

            if start >= min_start and end <= max_end:
                result = (
                    start,
                    end,
                )

                if with_reverse_map:
                    result.append(
                        list(
                            set(
                                [
                                    item
                                    for sublist in revmap[min(idx) : max(idx) + 1]
                                    for item in sublist
                                ]
                            )
                        )
                    )

                final.append(result)

    return final


def find_union(*intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Union all intervals, input can be any number of interval lists.

    Args:
        *intervals (List[Tuple[int, int]]): Lists containing intervals.

    Returns:
        List[Tuple[int, int]]: Union of all intervals.
    """
    all_intervals = []
    for intval in intervals:
        all_intervals.extend(intval)

    return find_unary_union(all_intervals, with_reverse_map=False)


def find_intersect(*intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Intersection of all intervals, input can be any number of interval lists.

    Args:
        *intervals (List[Tuple[int, int]]): List containing intervals.

    Returns:
        List[Tuple[int, int]]: Intersection across all intervals.
    """
    all_intervals = []
    for intval in intervals:
        all_intervals.extend(intval)

    return find_unary_intersect(
        all_intervals, with_reverse_map=False, threshold=len(intervals) - 1
    )


def find_diff(
    a: List[Tuple[int, int]], b: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """Set difference between intervals in a and b.

    Args:
        a, b (List[Tuple[int, int]]): Lists containing intervals.

    Returns:
        List[Tuple[int, int]]: Intersection of all intervals.
    """

    a_max_end = max([x[1] for x in a])
    b_max_end = max([x[1] for x in b])
    max_end = max(a_max_end, b_max_end)

    min_start = min([x[0] for x in a])

    a_intvals, _ = create_np_interval_vector(
        intervals=a, with_reverse_map=False, force_size=max_end, dont_sum=True
    )

    b_intvals, _ = create_np_interval_vector(
        intervals=b, with_reverse_map=False, force_size=max_end, dont_sum=True
    )

    diff_intvals = a_intvals[:a_max_end] - b_intvals[:a_max_end]

    groups = [
        (x[0], x[1])
        for x in zip(
            split(range(len(diff_intvals)), where(diff_intvals == 0)[0] + 1),
            split(diff_intvals, where(diff_intvals == 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if len(cov) > 0 and cov[0] > 0:
            start = min(idx) + 1
            end = max(idx) + 1 if (max(idx) == len(diff_intvals) - 1) else max(idx)

            if start >= min_start and end <= a_max_end:
                result = (start, end)
                final.append(result)

    return final


def compute_mean(
    intervals: List[Tuple[int, int]], values: Union[List[int], List[float]]
) -> Tuple[ndarray, ndarray]:
    """Compute coverage and sum ndarrays.

    Args:
        intervals (List[Tuple[int, int]]): Input interval list.
        values (Union[List[int], List[float]]): An ``int`` or ``float`` vector with
            the same size as intervals.

    Returns:
        Tuple[ndarray, ndarray]: A tuple with coverage and sum ndarrays.
    """
    max_end = None
    if max_end is None:
        max_end = max([x[1] for x in intervals])

    np_cov = zeros(max_end)
    np_sum = zeros(max_end)

    for idx in range(len(intervals)):
        i = intervals[idx]
        np_cov[i[0] - 1 : i[1]] += 1
        np_sum[i[0] - 1 : i[1]] += values[idx]

    return (np_cov, np_sum)


OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]


def find_overlaps(
    subject: List[Tuple[int, int]],
    query: List[Tuple[int, int]],
    max_gap: int = -1,
    min_overlap: int = 1,
    query_type: str = "any",
) -> List[Tuple[Tuple[int, int], int, List[int]]]:
    """Find overlaps between subject and query.

    Args:
        subject (List[Tuple[int, int]]): Intervals.
        query (List[Tuple[int, int]]): Query intervals.
        max_gap (int, optional): Maximum gap allowed. Defaults to -1 for no gaps.
        min_overlap (int, optional): Minimum overlap needed for intervals to be
            considered as overlapping. Defaults to 1.
        query_type (str, optional): Query type, one of any
            - "any": any overlap is good
            - "start": overlap at the beginning of the intervals
            - "end": must overlap at the end of the intervals
            - "within": Fully contain the query interval.

            Defaults to "any".

    Raises:
        ValueError: Query type is incorrect.

    Returns:
        List[Tuple[Tuple[int, int], int, List[int]]]: List of query intervals
        with their overlaps.
    """
    if query_type not in OVERLAP_QUERY_TYPES:
        raise ValueError(f"{query_type} must be one of {OVERLAP_QUERY_TYPES}")

    subject_ints, subject_revmap = create_np_interval_vector(
        intervals=subject, with_reverse_map=True
    )

    hits = []
    for idx in range(len(query)):
        q = query[idx]
        vec = subject_ints[q[0] - 1 + max_gap : q[1] + abs(max_gap)]
        overlaps = where(vec > 0)[0]

        if len(overlaps) < min_overlap:
            continue

        indices = subject_revmap[q[0] - 1 + max_gap : q[1] + abs(max_gap)]

        if query_type == "start":
            indices = indices[0]
        elif query_type == "end":
            indices = indices[-1]
        else:
            indices = list(set([item for sublist in indices for item in sublist]))
            if query_type == "within":
                tmp_idx = []
                for qidx in indices:
                    if subject[idx][0] <= q[0] and subject[idx][1] >= q[1]:
                        tmp_idx.append(qidx)

                indices = tmp_idx

        hits.append(((q[0], q[1]), idx + 1, indices))
    return hits


def find_nearest(
    subject: List[Tuple[int, int]],
    query: List[Tuple[int, int]],
    stepstart: int = 3,
    stepend: int = 3,
) -> List[Tuple[Tuple[int, int], int, List[int]]]:
    """Find nearest intervals in subject for even interval in query.

    Args:
        subject (List[Tuple[int, int]]): Intervals.
        query (List[Tuple[int, int]]): Query intervals.
        max_gap (int, optional): Maximum gap allowed. Defaults to -1 for no gaps.
        stepstart (int, optional): Step start. Defaults to 3.
        stepend (int, optional): Step end. Defaults to 3.

    Returns:
        List[Tuple[Tuple[int, int], int, List[int]]]: List of query intervals
        with their overlaps.
    """
    _, subject_revmap = create_np_interval_vector(
        intervals=subject, with_reverse_map=True
    )

    hits = []
    for idx in range(len(query)):
        q = query[idx]

        matches = 0
        counter = 0
        indices = []
        while matches == 0:
            slices = subject_revmap[
                q[0] - 1 - (counter * stepstart) : q[1] + (counter * stepend)
            ]

            indices = list(set([item for sublist in slices for item in sublist]))

            matches = len(indices)
            counter += 1

            if q[0] - 1 - (counter * stepstart) < 1 and q[1] + (
                counter * stepend
            ) > len(subject_revmap):
                matches = -1

            if stepend == 0 and q[0] - 1 - (counter * stepstart) < 1:
                matches = -1

            if stepstart == 0 and (counter * stepend) > len(subject_revmap):
                matches = -1

        hits.append(
            (
                (q[0], q[1]),
                idx + 1,
                indices,
                min(counter * stepstart, counter * stepend),
            )
        )

    return hits


def split_intervals(
    chrom: str, strand: str, start: int, end: int, step: int
) -> List[Tuple]:
    """Split an interval range into equal bins. pretty much a fancy range function. realizes the range.

    Args:
        chrom (str): chromosome name.
        strand (str): Strand information.
        start (int): Start interval.
        end (int): End interval.
        step (int): Width or step of each interval.

    Returns:
        List[Tuple]: List of internals split into bins.
    """
    bins = []
    for i in range(start, end + 1, step):
        bins.append((chrom, strand, i, min(i + step - 1, end)))

    return bins


def slide_intervals(
    chrom: str, strand: str, start: int, end: int, width: int, step: int
):
    """Sliding intervals. pretty much a fancy range function. realizes the range.

    Args:
        chrom (str): chromosome name.
        strand (str): Strand information.
        start (int): Start interval.
        end (int): End interval.
        step (int): Width or step of each interval.

    Returns:
        List[Tuple]: List of internals split into bins.
    """
    bins = []

    if end - width < start:
        bins.append((chrom, strand, start, end))
    else:
        for i in range(start, end - width + 2, step):
            bins.append((chrom, strand, i, min(i + width - 1, end)))

    return bins


def adjust_interval(
    row: MutableMapping[str, Any],
    shift_start: int = 0,
    shift_end: int = 0,
) -> Tuple[int, int]:
    """Shift genomic intervals for a ``row`` by their ``shift_start`` or ``shift_end``.

    Note: These values can be negative as well.

    ``start`` and ``end`` cannot be negative after shift! If so, they are set to 0.

    Args:
        row (MutableMapping[str, Any]): A row from `GenomicRanges`.
        shift_start (int, optional): Number of positions to shift start by.
            Defaults to 0.
        shift_end (int, optional): Number of positions to shift end by.
            Defaults to 0.

    Returns:
        Tuple[int, int]: Adjusted start and end interval.
    """
    return (max(0, row["starts"] + shift_start), max(0, row["ends"] + shift_end))
