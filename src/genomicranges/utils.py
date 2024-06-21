from typing import List, Optional, Sequence, Tuple, Union

import biocutils as ut
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

STRAND_MAP = {"+": 1, "-": -1, "*": 0}
REV_STRAND_MAP = {"1": "+", "-1": "-", "0": "*"}


def sanitize_strand_vector(
    strand: Union[Sequence[str], Sequence[int], np.ndarray]
) -> np.ndarray:
    """Create a numpy representation for ``strand``.

    Mapping: 1 for "+" (forward strand), 0 for "*" (any strand) and -1 for "-" (reverse strand).

    Args:
        strand: List of strand.

    Raises:
        ValueError:
            If strand is None.
            If strand contains values other than +,- and *.
            If strand is not a numpy vector, string of integers or strings.

    Returns:
        A numpy vector.
    """
    if strand is None:
        raise ValueError("'strand' cannot be None.")

    if isinstance(strand, np.ndarray):
        if len(strand.shape) > 1:
            raise ValueError("'strand' must be a 1-dimensional vector.")

        if not set(np.unique(strand)).issubset([-1, 0, 1]):
            raise ValueError(
                "'strand' must only contain values 1 (forward strand), -1 (reverse strand) or 0 (reverse strand)."
            )
        return strand.astype(np.int8)
    elif ut.is_list_of_type(strand, str):
        if not set(strand).issubset(["+", "-", "*"]):
            raise ValueError("Values in 'strand' must be either +, - or *.")
        return np.asarray([STRAND_MAP[x] for x in strand], dtype=np.int8)
    elif ut.is_list_of_type(strand, (int, float)):
        if not set(strand).issubset([1, 0, -1]):
            raise ValueError(
                "'strand' must only contain values 1 (forward strand), -1 (reverse strand) or 0 (reverse strand)."
            )
        return np.asarray(strand, dtype=np.int8)
    else:
        raise TypeError(
            "'strand' must be either a numpy vector, a list of integers or strings representing strand."
        )


def _sanitize_vec(x: Sequence):
    if isinstance(x, np.ma.MaskedArray):
        x.filled(fill_value=None)
        return x.tolist()

    return list(x)


def _sanitize_strand_search_ops(query_strand, subject_strand):
    query_strand = REV_STRAND_MAP[query_strand]
    subject_strand = REV_STRAND_MAP[subject_strand]

    out = None

    if query_strand == "+":
        if subject_strand == "+":
            out = "+"
        elif subject_strand == "-":
            out = None
        elif subject_strand == "*":
            out = "+"
    elif query_strand == "-":
        if subject_strand == "+":
            out = None
        elif subject_strand == "-":
            out = "-"
        elif subject_strand == "*":
            out = "-"
    elif query_strand == "*":
        if subject_strand == "*":
            out = "+"
        elif subject_strand == "-":
            out = "-"
        elif subject_strand == "*":
            out = "-"

    if out is None:
        return None

    return STRAND_MAP[out]


def split_intervals(start: int, end: int, step: int) -> List:
    """Split an interval range into equal bins.

    Args:
        start:
            Start interval.

        end:
            End interval.

        step:
            Width or step of each interval.

    Returns:
        List of intervals split into bins.
    """
    bins = []
    for i in range(start, end + 1, step):
        bins.append((i, min(i + step - 1, end) - i))

    return bins


def slide_intervals(start: int, end: int, width: int, step: int) -> List:
    """Sliding intervals.

    Args:
        start:
            Start interval.

        end:
            End interval.

        step:
            Step of each interval.

        width:
            Width of each interval.

    Returns:
        List of intervals split into bins.
    """
    bins = []

    if end - width < start:
        bins.append((start, end - start))
    else:
        for i in range(start, end - width + 2, step):
            bins.append((i, min(i + width - 1, end) - i))

    return bins


def create_np_vector(
    intervals: List[Tuple[int, int]],
    with_reverse_map: bool = False,
    force_size: Optional[int] = None,
    dont_sum: bool = False,
    value: int = 1,
) -> Tuple[np.ndarray, Optional[List]]:
    """Represent intervals and calculate coverage.

    Args:
        intervals:
            Input interval vector.

        with_reverse_map:
            Return map of indices? Defaults to False.

        force_size:
            Force size of the array.

        dont_sum:
            Do not sum. Defaults to False.

        value:
            Default value to increment. Defaults to 1.

    Returns:
        A numpy array representing coverage from the
        intervals and optionally a reverse index map.
    """
    if len(intervals) < 1:
        return intervals

    max_end = force_size
    if max_end is None:
        max_end = max([x[1] for x in intervals])
    cov = np.zeros(max_end)

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
