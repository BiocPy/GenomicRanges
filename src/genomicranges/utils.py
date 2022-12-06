from typing import MutableMapping, Any, Optional, Tuple


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


#
#  The following methods are for computing gaps between genomic regions
#
def calc_start_gap(
    row: MutableMapping[str, Any], name: Tuple[str, str], start_limit: int
) -> Tuple:
    """Give a genomic position and chromosome limits, calculate gap

    Args:
        row (MutableMapping[str, Any]): a row from GenomicRanges
        name (Tuple[str, str]): Tuple with chromosome name and strand
        start_limit (int): set the chromosome start

    Returns:
        Tuple: a tuple with the gap interval
    """
    interval = None

    if row["starts"] < start_limit:
        return interval

    interval = (name[0], name[1], start_limit, row["starts"] - 1)

    return interval


def calc_end_gap(
    row: MutableMapping[str, Any],
    name: Tuple[str, str],
    end_limit: Optional[int] = None,
) -> Tuple:
    """Give a genomic position and chromosome limits, calculate gap

    Args:
        row (MutableMapping[str, Any]): a row from GenomicRanges
        name (Tuple[str, str]): Tuple with chromosome name and strand
        end_limit (int): set the chromosome end 

    Returns:
        Tuple: a tuple with the gap interval
    """
    interval = None

    if end_limit is None or (end_limit is not None and row["ends"] > end_limit):
        return interval

    interval = (name[0], name[1], row["ends"] + 1, end_limit)

    return interval


def calc_between_gap(
    rowf: MutableMapping[str, Any],
    rowl: MutableMapping[str, Any],
    name: Tuple[str, str],
    start_limit: int,
    end_limit: Optional[int] = None,
) -> Tuple:
    """Give two genomic positions and chromosome limits, calculate gap

    Args:
        rowf (MutableMapping[str, Any]): a row from GenomicRanges
        rowl (MutableMapping[str, Any]): a row from GenomicRanges
        name (Tuple[str, str]): Tuple with chromosome name and strand
        start_limit (int): set the chromosome start
        end_limit (int): set the chromosome end 

    Returns:
        Tuple: a tuple with the gap interval
    """
    interval = None

    # just in case they are not sorted
    sorted_first = rowf
    sorted_last = rowl
    if rowf["starts"] > rowl["starts"]:
        sorted_first = rowl
        sorted_last = rowf

    # overlap condition
    if (
        sorted_first["starts"] < sorted_last["ends"]
        and sorted_first["ends"] > sorted_last["starts"]
    ):
        return interval

    gap_start = sorted_first["ends"] + 1
    if gap_start < start_limit:
        return interval

    gap_end = sorted_last["starts"] - 1
    if end_limit is not None and gap_end > end_limit:
        return interval

    interval = (name[0], name[1], gap_start, gap_end)

    return interval
