from iranges import IRanges

from genomicranges import GenomicRanges

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

gr = GenomicRanges(
    seqnames=[
        "chr1",
        "chr2",
        "chr3",
        "chr2",
        "chr3",
    ],
    ranges=IRanges([x for x in range(101, 106)], [11, 21, 25, 30, 5]),
    strand=["*", "-", "*", "+", "-"],
)


def test_tile():
    assert gr is not None

    tiles = gr.tile(n=2)
    assert tiles is not None
    assert isinstance(tiles, list)
    assert len(tiles) == 5
    assert sum(len(x) for x in tiles) == 10

    assert tiles[0].get_seqnames() == ["chr1", "chr1"]


def test_slide():
    assert gr is not None

    tiles = gr.sliding_windows(width=3)

    assert tiles is not None
    assert isinstance(tiles, list)
    assert sum(len(x) for x in tiles) == 82


def test_tile_genome():
    seqlengths = {"chr1": 60, "chr2": 20, "chr3": 25}

    tiles = GenomicRanges.tile_genome(seqlengths=seqlengths, ntile=5)

    assert tiles is not None
    assert isinstance(tiles, GenomicRanges)
    assert len(tiles) == 7
