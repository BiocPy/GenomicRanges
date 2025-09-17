from itertools import groupby
from typing import List, Sequence, Union

import biocutils as ut
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

STRAND_MAP = {"+": 1, "-": -1, "*": 0}
REV_STRAND_MAP = {"1": "+", "-1": "-", "0": "*"}


def sanitize_strand_vector(strand: Union[Sequence[str], Sequence[int], np.ndarray]) -> np.ndarray:
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
        raise TypeError("'strand' must be either a numpy vector, a list of integers or strings representing strand.")


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


def group_by_indices(groups: list) -> dict:
    return {k: [x[0] for x in v] for k, v in groupby(sorted(enumerate(groups), key=lambda x: x[1]), lambda x: x[1])}


def compute_up_down(starts, ends, strands, upstream, downstream, site: str = "TSS"):
    """Compute promoter or terminator regions for genomic ranges.

    Args:
        x:
            GenomicRanges object

        upstream:
            Number of bases upstream (scalar or array)

        downstream:
            Number of bases downstream (scalar or array)

        site:
            "TSS" for promoters or "TES" for terminators

    Returns:
        New starts and ends.
    """
    if isinstance(upstream, (int, float)):
        upstream = np.full(len(starts), upstream, dtype=np.int32)
    if isinstance(downstream, (int, float)):
        downstream = np.full(len(starts), downstream, dtype=np.int32)

    minus_mask = strands == -1

    new_starts = np.zeros_like(starts)
    new_ends = np.zeros_like(ends)

    if site == "TSS":
        # For "+" or "*" strand
        plus_mask = ~minus_mask
        new_starts[plus_mask] = starts[plus_mask] - upstream[plus_mask]
        new_ends[plus_mask] = starts[plus_mask] + downstream[plus_mask] - 1

        # For "-" strand
        new_starts[minus_mask] = ends[minus_mask] - downstream[minus_mask] + 1
        new_ends[minus_mask] = ends[minus_mask] + upstream[minus_mask]
    else:  # TES
        # For "+" or "*" strand
        plus_mask = ~minus_mask
        new_starts[plus_mask] = ends[plus_mask] - upstream[plus_mask]
        new_ends[plus_mask] = ends[plus_mask] + downstream[plus_mask] - 1

        # For "-" strand
        new_starts[minus_mask] = starts[minus_mask] - downstream[minus_mask] + 1
        new_ends[minus_mask] = starts[minus_mask] + upstream[minus_mask]

    return new_starts, new_ends - new_starts + 1


def extract_groups_from_granges(x, ignore_strand=False):
    if ignore_strand:
        groups = []
        for idx, seq in enumerate(x._seqinfo._seqnames):
            matches = np.where(x._seqnames == idx)[0]
            groups.append((seq, matches))
        return groups
    else:
        ## TODO: needs some rethinking to speed this up
        combined = np.stack((x._seqnames, x._strand), axis=-1)
        unique_groups, inverse_indices = np.unique(combined, return_inverse=True, axis=0)

        groups = []
        for i, grp_array in enumerate(unique_groups):
            seq_idx, strand_val = grp_array[0], grp_array[1]
            seq_name = x._seqinfo._seqnames[seq_idx]

            matches = np.where(inverse_indices == i)[0]
            groups.append(((seq_name, strand_val), matches))

        return groups


def wrapper_follow_precede(args):
    """Processes a single group for precede/follow operations.
    This function is designed to be called by a multiprocessing pool.
    """
    (
        method_name,
        s_group,
        q_group,
        self_ranges,
        query_ranges,
        self_strand,
        query_strand,
        ignore_strand,
    ) = args

    method_to_call = getattr(self_ranges[s_group], method_name)
    res_idx = method_to_call(query=query_ranges[q_group], select="all")

    _q_hits = np.asarray([q_group[j] for j in res_idx.get_column("query_hits")])
    _s_hits = np.asarray([s_group[x] for x in res_idx.get_column("self_hits")])

    if not ignore_strand:
        s_strands = self_strand[s_group][res_idx.get_column("self_hits")]
        q_strands = query_strand[q_group][res_idx.get_column("query_hits")]

        mask = (s_strands == q_strands) | (s_strands == 0) | (q_strands == 0)

        _q_hits = _q_hits[mask]
        _s_hits = _s_hits[mask]

    return _q_hits, _s_hits
