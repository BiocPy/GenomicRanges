from typing import Sequence, Union

import biocutils as ut
import numpy as np

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

STRAND_MAP = {"+": 1, "-": -1, "*": 0}


def make_strand_vector(strand: Union[Sequence[str], Sequence[int]]) -> np.ndarray:
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

    if isinstance(strand, ut.StringList):
        if not set(strand).issubset(["+", "-", "+"]):
            raise ValueError("Values in 'strand' must be either +, - or *.")
        return np.ndarray([STRAND_MAP[x] for x in strand])
    elif ut.is_list_of_type(strand, int):
        return np.ndarray(strand)
    else:
        TypeError(
            "'strand' must be either a numpy vector, a list of integers or strings representing strand."
        )
