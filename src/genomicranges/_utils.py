# from typing import MutableMapping, Any, Optional, Tuple, List, Union
# from functools import cmp_to_key
# import numpy as np
# from random import random

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

# # credits: with help from SO - https://stackoverflow.com/a/31601579
# def calc_disjoint_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
#     """split intervals into a set of disjoint non-overlapping intervals

#     Args:
#         intervals (List[Tuple[int, int]]): List of input intervals.

#     Returns:
#         List[Tuple[int, int]]: List of disjoint intervals from input.
#     """
#     if len(intervals) < 1:
#         return intervals

#     disjoint = []

#     # get boundaries
#     all_starts = [(x[0], 1) for x in intervals]
#     all_ends = [(x[1] + 1, -1) for x in intervals]

#     delta = None
#     split_counter = 0

#     all_points = all_starts + all_ends

#     def compare(x, y):
#         return x[0] - y[0]

#     all_points = sorted(all_points, key=cmp_to_key(compare))

#     disjoint = []
#     for tint in all_points:
#         if delta is not None and tint[0] > delta and split_counter != 0:
#             disjoint.append((delta, tint[0] - 1))

#         delta = tint[0]
#         split_counter += tint[1]

#     return disjoint
