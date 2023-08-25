from .GenomicRanges import GenomicRanges
from collections import UserDict

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class GenomicRangesList(UserDict):
    def __setitem__(self, key, value):
        if not isinstance(value, GenomicRanges):
            raise TypeError("Value must be of type `GenomicRanges`.")

        super().__setitem__(key, value)
