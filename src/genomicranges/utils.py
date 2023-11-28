from typing import List, Sequence, Union

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
        return strand

    if ut.is_list_of_type(strand, str):
        if not set(strand).issubset(["+", "-", "*"]):
            raise ValueError("Values in 'strand' must be either +, - or *.")
        return np.array([STRAND_MAP[x] for x in strand])
    elif ut.is_list_of_type(strand, int):
        return np.array(strand)
    else:
        TypeError(
            "'strand' must be either a numpy vector, a list of integers or strings representing strand."
        )


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
