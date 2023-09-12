from typing import Optional, Sequence, Tuple, Union

from biocframe import BiocFrame
from pandas import DataFrame

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

SlicerTypes = Union[Sequence[int], Sequence[bool], Sequence[str], slice, int, str]
SlicerArgTypes = Union[Sequence[str], Tuple[SlicerTypes, Optional[SlicerTypes]]]
BiocOrPandasFrame = Union[DataFrame, BiocFrame]
