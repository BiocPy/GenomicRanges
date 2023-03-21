from genomicranges.utils import (
    find_union,
    find_gaps,
    find_disjoin,
    find_intersect,
    find_diff,
)

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


a = [(101, 112), (113, 128)]
b = [(101, 111), (105, 125), (136, 140)]

# intervals = [(1, 1), (1, 4), (2, 2), (4, 4), (5, 5), (6, 8), (7, 9), (10, 10)]


def test_union():
    res = find_union(a, b)

    assert res is not None
    assert len(res) == 2
    assert res == [(101, 128), (136, 140)]


def test_gaps():
    intervals = a + b
    res = find_gaps(intervals, start_limit=None)

    assert res is not None
    assert len(res) == 1
    assert res == [(129, 135)]


def test_disjoin():

    intervals = a + b
    res1 = find_disjoin(intervals)
    # res2 = calc_disjoint_intervals(intervals)

    # assert len(res2) == 6
    assert len(res1) >= 5

    # assert res1 == res2


def test_intersect():
    res = find_intersect(a, b)

    assert res is not None
    assert len(res) == 1
    assert res == [(101, 128)]


def test_diff():
    res = find_diff(a, b)

    assert res is not None
    assert len(res) == 1
    assert res == [(126, 128)]
