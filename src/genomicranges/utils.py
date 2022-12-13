from typing import MutableMapping, Any, Optional, Tuple, List, Union
from functools import cmp_to_key
import numpy as np
from random import random


def calc_row_gapwidth(
    a: MutableMapping[str, Any], b: MutableMapping[str, Any]
) -> Optional[int]:
    """Given two genomic positions from GenomicRanges, calculates gapwidth
        returns 
        - `nan` if the sequences are not comparable
        - 0 if they overlap
        - a number if there is a gap

    Args:
        a (MutableMapping[str, Any]): a row from GenomicRanges
        b (MutableMapping[str, Any]): a row from GenomicRanges

    Returns:
        Optional[int]: gapwidth if exists else None
    """
    if a["seqnames"] != b["seqnames"] or a["strand"] != b["strand"] != "-":
        return float("nan")

    if (a["starts"] <= b["starts"] <= a["ends"]) or (
        a["starts"] <= b["ends"] <= a["ends"]
    ):
        return 0

    return b["starts"] - a["ends"]


def np_split_groups_ne(
    ary: np.ndarray, step: int = 0
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Split a numpy array by consecutive values

    Args:
        ary (np.ndarray): a numpy array
        step (int, Optional): value to split consecutive intervals by

    Returns:
        Tuple(List[np.ndarray], List[np.ndarray]): tuple with indices and their coverage
    """
    return [
        (x[0], x[1])
        for x in zip(
            np.split(range(len(ary)), np.where(np.diff(ary) != step)[0] + 1),
            np.split(ary, np.where(np.diff(ary) != step)[0] + 1),
        )
    ]


def create_np_interval_vector(
    intervals: List[Tuple[int, int]],
    withRevMap: bool = False,
    forceSize: Optional[int] = None,
    dontSum: bool = False,
    value: int = 1,
) -> Tuple[np.ndarray, Optional[List]]:
    """Represent intervals/calculate coverage

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices. Defaults to False.
        forceSize (Optional[int], optional): force size of the array.
        dontSum (bool, optional): do not sum. Defaults to False.
        value (int, optional): default value to increment. Defaults to 1.

    Returns:
        Tuple[np.ndarray, Optional[List]]: a numpy array representing coverage from the intervals
    """
    if len(intervals) < 1:
        return intervals

    max_end = forceSize
    if max_end is None:
        max_end = max([x[1] for x in intervals])
    cov = np.zeros(max_end)

    revmap = None
    if withRevMap:
        revmap = [[] for _ in range(max_end)]

    for idx in range(len(intervals)):
        i = intervals[idx]

        if dontSum:
            cov[i[0] - 1 : i[1]] = value
        else:
            cov[i[0] - 1 : i[1]] += value

        if withRevMap:
            _ = [revmap[x].append(idx + 1) for x in range(i[0] - 1, i[1])]

    return cov, revmap


def find_unary_union(
    intervals: List[Tuple[int, int]], withRevMap: bool = False
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Union of all overlapping intervals

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices. Defaults to False.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List containing tuples with 
            (start, end, indices) for each contiguous region
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, withRevMap=withRevMap
    )

    # splitting intervals by where gaps=0 gives us the union of all intervals
    groups = [
        (x[0], x[1])
        for x in zip(
            np.split(range(len(np_intvals)), np.where(np_intvals == 0)[0] + 1),
            np.split(np_intvals, np.where(np_intvals == 0)[0] + 1),
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

            if withRevMap:
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
    """Find Gaps across all intervals

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        start_limit (int, optional): start for identifying gaps. Defaults to 1.
        end_limit (int, optional): end interval. Defaults to None.

    Returns:
        List[Tuple[int, int]]: List of tuples with 
            (start, end) for each gap
    """

    min_start = start_limit
    if min_start is None:
        min_start = min([x[0] for x in intervals])

    max_end = end_limit
    if max_end is None:
        max_end = max(x[1] for x in intervals)

    np_intvals, _ = create_np_interval_vector(
        intervals=intervals, withRevMap=False, forceSize=max_end
    )

    # splitting intervals by diff gives us the gaps
    groups = [
        (x[0], x[1])
        for x in zip(
            np.split(range(len(np_intvals)), np.where(np.diff(np_intvals) != 0)[0] + 1),
            np.split(np_intvals, np.where(np.diff(np_intvals) != 0)[0] + 1),
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
    intervals: List[Tuple[int, int]], withRevMap: bool = False
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Find disjoin sets across all intervals

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices. Defaults to False.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List of tuples with 
            (start, end)
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, withRevMap=withRevMap
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
            np.split(range(len(np_intvals)), np.where(np.diff(np_intvals) != 0)[0] + 1),
            np.split(np_intvals, np.where(np.diff(np_intvals) != 0)[0] + 1),
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

                if withRevMap:
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
    intervals: List[Tuple[int, int]], withRevMap: bool = False, threshold: int = 1
) -> Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]:
    """Intersection of all intervals

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices?. Defaults to False.
        threshold (int, optional): threshold for cutoff. Defaults to 1.

    Returns:
        Union[List[Tuple[int, int, Optional[List[int]]]], List[Tuple[int, int]]]: List containing tuples with 
            (start, end, indices) for each contiguous region
    """
    np_intvals, revmap = create_np_interval_vector(
        intervals=intervals, withRevMap=withRevMap
    )

    min_start = min([x[0] for x in intervals])
    max_end = max(x[1] for x in intervals)

    # splitting intervals by where gaps=0 gives us the union of all intervals
    groups = [
        (x[0], x[1])
        for x in zip(
            np.split(range(len(np_intvals)), np.where(np_intvals == 0)[0] + 1),
            np.split(np_intvals, np.where(np_intvals == 0)[0] + 1),
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

                if withRevMap:
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
    """Union all intervals, input can be any number of interval lists

    Args:
        *intervals (List[Tuple[int, int]]): lists containing intervals

    Returns:
        List[Tuple[int, int]]: Union of all intervals
    """
    all_intervals = []
    for intval in intervals:
        all_intervals.extend(intval)

    return find_unary_union(all_intervals, withRevMap=False)


def find_intersect(*intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Intersection of all intervals, input can be any number of interval lists

    Args:
        *intervals (List[Tuple[int, int]]): lists containing intervals

    Returns:
        List[Tuple[int, int]]: Intersection of all intervals
    """
    all_intervals = []
    for intval in intervals:
        all_intervals.extend(intval)

    return find_unary_intersect(
        all_intervals, withRevMap=False, threshold=len(intervals) - 1
    )


def find_diff(
    a: List[Tuple[int, int]], b: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """Set difference between intervals in a and b

    Args:
        a, b (List[Tuple[int, int]]): lists containing intervals

    Returns:
        List[Tuple[int, int]]: Intersection of all intervals
    """

    a_max_end = max([x[1] for x in a])
    b_max_end = max([x[1] for x in b])
    max_end = max(a_max_end, b_max_end)

    min_start = min([x[0] for x in a])

    a_intvals, _ = create_np_interval_vector(
        intervals=a, withRevMap=False, forceSize=max_end, dontSum=True
    )

    b_intvals, _ = create_np_interval_vector(
        intervals=b, withRevMap=False, forceSize=max_end, dontSum=True
    )

    diff_intvals = a_intvals[:a_max_end] - b_intvals[:a_max_end]

    groups = [
        (x[0], x[1])
        for x in zip(
            np.split(range(len(diff_intvals)), np.where(diff_intvals == 0)[0] + 1),
            np.split(diff_intvals, np.where(diff_intvals == 0)[0] + 1),
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
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute coverage and sum ndarrays

    Args:
        intervals (List[Tuple[int, int]]): input interval list
        values (Union[List[int], List[float]]): an int or float vector with the same size as intervals

    Returns:
        Tuple[np.ndarray, np.ndarray]: a tuple with coverage and sum noarrays
    """
    max_end = None
    if max_end is None:
        max_end = max([x[1] for x in intervals])

    np_cov = np.zeros(max_end)
    np_sum = np.zeros(max_end)

    for idx in range(len(intervals)):
        i = intervals[idx]
        np_cov[i[0] - 1 : i[1]] += 1
        np_sum[i[0] - 1 : i[1]] += values[idx]

    return (np_cov, np_sum)


def compute_mean(
    intervals: List[Tuple[int, int]], values: Union[List[int], List[float]]
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute coverage and sum ndarrays

    Args:
        intervals (List[Tuple[int, int]]): input interval list
        values (Union[List[int], List[float]]): an int or float vector with the same size as intervals

    Returns:
        Tuple[np.ndarray, np.ndarray]: a tuple with coverage and sum noarrays
    """
    max_end = None
    if max_end is None:
        max_end = max([x[1] for x in intervals])

    np_cov = np.zeros(max_end)
    np_sum = np.zeros(max_end)

    for idx in range(len(intervals)):
        i = intervals[idx]
        np_cov[i[0] - 1 : i[1]] += 1
        np_sum[i[0] - 1 : i[1]] += values[idx]

    return (np_cov, np_sum)


OVERLAP_QUERY_TYPES = ["any", "start", "end", "within"]


def find_overlaps(
    subject: List[Tuple[int, int]],
    query: List[Tuple[int, int]],
    maxGap: int = -1,
    minOverlap: int = 1,
    queryType: str = "any",
) -> List[Tuple[Tuple[int, int], int, List[int]]]:
    """Find overlaps between subject and query

    Args:
        subject (List[Tuple[int, int]]): intervals
        query (List[Tuple[int, int]]): query intervals
        maxGap (int, optional): maximum gap allowed. Defaults to -1 for no gaps.
        minOverlap (int, optional): minimum overlap needed for intervals to be considered as overlapping. Defaults to 1.
        queryType (str, optional): query type, one of any
                "any": any overlap is good
                "start": overlap at the beginning of the intervals
                "end": must overlap at the end of the intervals
                "within": Fully contain the query interval. 
            Defaults to "any".

    Raises:
        ValueError: query type is incorrect

    Returns:
        List[Tuple[Tuple[int, int], int, List[int]]]: list of query intervals with their overlaps
    """
    if queryType not in OVERLAP_QUERY_TYPES:
        raise ValueError(f"{queryType} must be one of {OVERLAP_QUERY_TYPES}")

    subject_ints, subject_revmap = create_np_interval_vector(
        intervals=subject, withRevMap=True
    )

    hits = []
    for idx in range(len(query)):
        q = query[idx]
        vec = subject_ints[q[0] - 1 + maxGap : q[1] + abs(maxGap)]
        overlaps = np.where(vec > 0)[0]

        if len(overlaps) < minOverlap:
            continue

        indices = subject_revmap[q[0] - 1 + maxGap : q[1] + abs(maxGap)]

        if queryType == "start":
            indices = indices[0]
        elif queryType == "end":
            indices = indices[-1]
        else:
            indices = list(set([item for sublist in indices for item in sublist]))
            if queryType == "within":
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
    """Find nearest intervals in subject for even interval in query

    Args:
        subject (List[Tuple[int, int]]): intervals
        query (List[Tuple[int, int]]): query intervals
        maxGap (int, optional): maximum gap allowed. Defaults to -1 for no gaps.
        stepstart (int, optional): step start. Defaults to 3.
        stepend (int, optional): step end. Defaults to 3.

    Returns:
        List[Tuple[Tuple[int, int], int, List[int]]]: list of query intervals with their overlaps
    """
    _, subject_revmap = create_np_interval_vector(intervals=subject, withRevMap=True)

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


def split_intervals(chrom: str, strand: str, start: int, end: int, step: int):
    """split an interval range into equal bins. pretty much a fancy range function.
        realizes the range.

    Args:
        chrom (str): chrom
        strand (str): strand
        start (int): start interval
        end (int): end interval
        step (int): width or step of each interval
    """
    bins = []
    for i in range(start, end + 1, step):
        bins.append((chrom, strand, i, min(i + step - 1, end)))

    return bins


def slide_intervals(
    chrom: str, strand: str, start: int, end: int, width: int, step: int
):
    """Sliding intervals. pretty much a fancy range function.
        realizes the range.

    Args:
        chrom (str): chrom
        strand (str): strand
        start (int): start interval
        end (int): end interval
        step (int): width or step of each interval
    """
    bins = []

    if end - width < start:
        bins.append((chrom, strand, start, end))
    else:
        for i in range(start, end - width + 2, step):
            bins.append((chrom, strand, i, min(i + width - 1, end)))

    return bins
