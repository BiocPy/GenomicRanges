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


# #
# #  The following methods are for computing gaps between genomic regions
# #
# def calc_start_gap(
#     row: MutableMapping[str, Any], name: Tuple[str, str], start_limit: int
# ) -> Tuple:
#     """Give a genomic position and chromosome limits, calculate gap

#     Args:
#         row (MutableMapping[str, Any]): a row from GenomicRanges
#         name (Tuple[str, str]): Tuple with chromosome name and strand
#         start_limit (int): set the chromosome start

#     Returns:
#         Tuple: a tuple with the gap interval
#     """
#     interval = None

#     if row["starts"] < start_limit:
#         return interval

#     interval = (name[0], name[1], start_limit, row["starts"] - 1)

#     return interval


# def calc_end_gap(
#     row: MutableMapping[str, Any],
#     name: Tuple[str, str],
#     end_limit: Optional[int] = None,
# ) -> Tuple:
#     """Give a genomic position and chromosome limits, calculate gap

#     Args:
#         row (MutableMapping[str, Any]): a row from GenomicRanges
#         name (Tuple[str, str]): Tuple with chromosome name and strand
#         end_limit (int): set the chromosome end 

#     Returns:
#         Tuple: a tuple with the gap interval
#     """
#     interval = None

#     if end_limit is None or (end_limit is not None and row["ends"] > end_limit):
#         return interval

#     interval = (name[0], name[1], row["ends"] + 1, end_limit)

#     return interval


# def calc_between_gap(
#     rowf: MutableMapping[str, Any],
#     rowl: MutableMapping[str, Any],
#     name: Tuple[str, str],
#     start_limit: int,
#     end_limit: Optional[int] = None,
# ) -> Tuple:
#     """Give two genomic positions and chromosome limits, calculate gap

#     Args:
#         rowf (MutableMapping[str, Any]): a row from GenomicRanges
#         rowl (MutableMapping[str, Any]): a row from GenomicRanges
#         name (Tuple[str, str]): Tuple with chromosome name and strand
#         start_limit (int): set the chromosome start
#         end_limit (int): set the chromosome end 

#     Returns:
#         Tuple: a tuple with the gap interval
#     """
#     interval = None

#     # just in case they are not sorted
#     sorted_first = rowf
#     sorted_last = rowl
#     if rowf["starts"] > rowl["starts"]:
#         sorted_first = rowl
#         sorted_last = rowf

#     # overlap condition
#     if (
#         sorted_first["starts"] < sorted_last["ends"]
#         and sorted_first["ends"] > sorted_last["starts"]
#     ):
#         return interval

#     gap_start = sorted_first["ends"] + 1
#     if gap_start < start_limit:
#         return interval

#     gap_end = sorted_last["starts"] - 1
#     if end_limit is not None and gap_end > end_limit:
#         return interval

#     interval = (name[0], name[1], gap_start, gap_end)

#     return interval


# def calc_gaps(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
#     """Calculate gaps in the intervals

#     Args:
#         intervals (List[Tuple[int, int]]): List of input intervals.

#     Returns:
#         List[Tuple[int, int]]: List of non-overlapping intervals/gaps from input.
#     """
#     if len(intervals) == 0:
#         return intervals

#     gaps = []
#     # just in case they are not already sorted, sort by start
#     intervals.sort(key=lambda x: x[0])

#     # Iterate over all the interval
#     for idx in range(1, len(intervals)):
#         prev_end = intervals[idx - 1][1]
#         curr_start = intervals[idx][0]

#         if prev_end < curr_start:
#             gaps.append([prev_end, curr_start])

#     return gaps


# def calc_disjoint_intervals(
#     intervals: List[Tuple[int, int]], indexMap: bool = False
# ) -> List[Tuple[int, int]]:
#     """Calculate disjoint intervals, uses the max disjoint algorithm

#     Args:
#         intervals (List[Tuple[int, int]]): List of input intervals.

#     Returns:
#         List[Tuple[int, int]]: List of disjoint intervals from input.
#     """

#     if len(intervals) < 1:
#         return intervals

#     # sort by ends
#     intervals.sort(key=lambda x: x[1])

#     disjoints = []
#     disjoints.append(intervals[0])

#     last_end = intervals[0][1]

#     for idx in range(1, len(intervals)):
#         current_start = intervals[idx][0]
#         current_end = intervals[idx][1]

#         if current_start > last_end:
#             intervals.append((current_start, current_end))
#             last_end = current_end

#     return disjoints


# credits: with help from SO - https://stackoverflow.com/a/31601579
def calc_disjoint_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """split intervals into a set of disjoint non-overlapping intervals

    Args:
        intervals (List[Tuple[int, int]]): List of input intervals.

    Returns:
        List[Tuple[int, int]]: List of disjoint intervals from input.
    """
    if len(intervals) < 1:
        return intervals

    disjoint = []

    # get boundaries
    all_starts = [(x[0], 1) for x in intervals]
    all_ends = [(x[1] + 1, -1) for x in intervals]

    delta = None
    split_counter = 0

    all_points = all_starts + all_ends

    def compare(x, y):
        return x[0] - y[0]

    all_points = sorted(all_points, key=cmp_to_key(compare))

    disjoint = []
    for tint in all_points:
        if delta is not None and tint[0] > delta and split_counter != 0:
            disjoint.append((delta, tint[0] - 1))

        delta = tint[0]
        split_counter += tint[1]

    return disjoint


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
) -> Tuple[np.ndarray, Optional[List]]:
    """Represent intervals as numpy vector

    Args:
        intervals (List[Tuple[int, int]]): input interval vector
        withRevMap (bool, optional): return map of indices. Defaults to False.
        forceSize (Optional[int], optional): force size of the array.

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

    print("start_limit, end_limit", start_limit, end_limit)
    min_start = start_limit
    if min_start is None:
        min_start = min([x[0] for x in intervals])

    max_end = end_limit
    if max_end is None:
        max_end = max(x[1] for x in intervals)

    np_intvals, _ = create_np_interval_vector(intervals=intervals, withRevMap=False, forceSize=max_end)

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

    a_max_end = max(x[1] for x in a)
    b_max_end = max(x[1] for x in b)
    max_end = max(a_max_end, b_max_end)

    a_intvals, _ = create_np_interval_vector(
        intervals=a, withRevMap=False, forceSize=max_end
    )

    b_intvals, _ = create_np_interval_vector(
        intervals=b, withRevMap=False, forceSize=max_end
    )

    diff_intvals = a_intvals[:a_max_end] - b_intvals[:a_max_end]

    groups = [
        (x[0], x[1])
        for x in zip(
            np.split(
                range(len(diff_intvals)), np.where(np.diff(diff_intvals) > 0)[0] + 1
            ),
            np.split(diff_intvals, np.where(np.diff(diff_intvals) > 0)[0] + 1),
        )
    ]

    final = []
    for idx, cov in groups:
        if cov[0] > 0:
            start = min(idx) + 1
            end = max(idx) + 1

            result = (start, end)
            final.append(result)

    return final
