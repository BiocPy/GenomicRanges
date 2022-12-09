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
) -> Tuple[np.ndarray, Optional[List]]:
    """Represent intervals as numpy vector

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices. Defaults to False.
        forceSize (Optional[int], optional): force size of the array.
        dontSum (bool, optional): do not sum

    Returns:
        Tuple[np.ndarray, Optional[List]]: a numpy array representing the intervals
    """
    if len(intervals) < 1:
        return intervals

    max_end = forceSize
    if max_end is None:
        max_end = max([x[1] for x in intervals])
    np_intvals = np.zeros(max_end)

    revmap = None
    if withRevMap:
        revmap = [[] for _ in range(max_end)]

    for idx in range(len(intervals)):
        i = intervals[idx]

        if dontSum:
            np_intvals[i[0] - 1 : i[1]] = 1
        else:
            np_intvals[i[0] - 1 : i[1]] += 1

        if withRevMap:
            tmp = [revmap[x].append(idx + 1) for x in range(i[0] - 1, i[1])]

    return np_intvals, revmap


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


def compute_mean(intervals: List[Tuple[int, int]], values: List[int]) -> np.ndarray:
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
