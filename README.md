[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)
[![PyPI-Server](https://img.shields.io/pypi/v/GenomicRanges.svg)](https://pypi.org/project/GenomicRanges/)
![Unit tests](https://github.com/BiocPy/GenomicRanges/actions/workflows/pypi-test.yml/badge.svg)

# GenomicRanges

GenomicRanges provides container classes designed to represent genomic locations and support genomic analysis. It is similar to Bioconductor's [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html). **_Intervals are inclusive on both ends and starts at 1._**

**Note: V0.4.0 is a complete overhaul of the package, as such the constructor to GenomicRanges has changed!**

To get started, install the package from [PyPI](https://pypi.org/project/genomicranges/)

```shell
pip install genomicranges
```

## `GenomicRanges`

`GenomicRanges` is the base class to represent and operate over genomic regions and annotations.

### From UCSC or GTF file

You can easily download and parse genome annotations from UCSC or load a genome annotation from a GTF file,

```python
import genomicranges

gr = genomicranges.read_gtf(<PATH TO GTF>)
# OR
gr = genomicranges.read_ucsc(genome="hg19")

print(gr)
## output
## GenomicRanges with 1760959 intervals & 10 metadata columns.
## ... truncating the console print ...
```

### from `IRanges` (Preferred way)

If you have all relevant information to create a GenomicRanges object

```python
from genomicranges import GenomicRanges
from iranges import IRanges
from biocframe import BiocFrame
from random import random

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
    mcols=BiocFrame(
        {
            "score": range(0, 5),
            "GC": [random() for _ in range(5)],
        }
    ),
)

print(gr)
```

    ## output
    GenomicRanges with 5 ranges and 5 metadata columns
              seqnames    ranges           strand    score                  GC
      <ndarray[int64]> <IRanges> <ndarray[int64]>   <list>              <list>
    0             chr1 101 - 112                * |      0  0.9577094406822213
    1             chr2 102 - 103                - |      1  0.8947394912063824
    2             chr1 103 - 128                * |      2  0.8635855064433551
    3             chr3 104 - 134                + |      3 0.20372412804486073
    4             chr2 109 - 111                - |      4  0.6892415071280414

### Pandas DataFrame

A common representation in Python is a pandas `DataFrame` for all tabular datasets. `DataFrame` must contain columns "seqnames", "starts", and "ends" to represent genomic intervals. Here's an example:

```python
import genomicranges import GenomicRanges
import pandas as pd

df = pd.DataFrame(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3", "chr2"],
        "starts": [101, 102, 103, 104, 109],
        "ends": [112, 103, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

gr = GenomicRanges.from_pandas(df)
print(gr)
```

    ## output
    GenomicRanges with 5 ranges and 5 metadata columns
              seqnames    ranges           strand    score                  GC
      <ndarray[int64]> <IRanges> <ndarray[int64]>   <list>              <list>
    0             chr1 101 - 112                * |      0  0.9577094406822213
    1             chr2 102 - 103                - |      1  0.8947394912063824
    2             chr1 103 - 128                * |      2  0.8635855064433551
    3             chr3 104 - 134                + |      3 0.20372412804486073
    4             chr2 109 - 111                - |      4  0.6892415071280414

### Interval Operations

`GenomicRanges` supports most [interval based operations](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).

```python
subject = genomicranges.from_ucsc(genome="hg38")

query = genomicranges.from_pandas(
    pd.DataFrame(
        {
            "seqnames": ["chr1", "chr2", "chr3"],
            "starts": [100, 115, 119],
            "ends": [103, 116, 120],
        }
    )
)

hits = subject.nearest(query)
print(hits)
```

    ## output
    BiocFrame with 3 rows and 2 columns
        query                subject_hits
        <list>                      <list>
    [0]      0                      [0, 1]
    [1]      1 [1677082, 1677083, 1677084]
    [2]      2          [1003411, 1003412]

## `GenomicRangesList`

Just as it sounds, a `GenomicRangesList` is a named-list like object. If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple genomic elements, usually where the genes start and end. Genes are themselves made of many sub-regions, e.g. exons. `GenomicRangesList` allows us to represent this nested structure.

**Currently, this class is limited in functionality.**

To construct a GenomicRangesList

```python
gr1 = GenomicRanges(
    seqnames=["chr1", "chr2", "chr1", "chr3"],
    ranges=IRanges([1, 3, 2, 4], [10, 30, 50, 60]),
    strand=["-", "+", "*", "+"],
    mcols=BiocFrame({"score": [1, 2, 3, 4]}),
)

gr2 = GenomicRanges(
    seqnames=["chr2", "chr4", "chr5"],
    ranges=IRanges([3, 6, 4], [30, 50, 60]),
    strand=["-", "+", "*"],
    mcols=BiocFrame({"score": [2, 3, 4]}),
)
grl = GenomicRangesList(ranges=[gr1, gr2], names=["gene1", "gene2"])
print(grl)
```

    ## output
    GenomicRangesList with 2 genomic elements

    Name: gene1
    GenomicRanges with 4 ranges and 4 metadata columns
                seqnames    ranges           strand    score
        <ndarray> <IRanges> <ndarray>   <list>
    [0]             chr1    1 - 11                - |      1
    [1]             chr2    3 - 33                + |      2
    [2]             chr1    2 - 52                * |      3
    [3]             chr3    4 - 64                + |      4
    ------
    seqinfo(3 sequences): chr1 chr2 chr3
    Name: gene2
    GenomicRanges with 3 ranges and 3 metadata columns
                seqnames    ranges           strand    score
        <ndarray> <IRanges> <ndarray>   <list>
    [0]             chr2    3 - 33                - |      2
    [1]             chr4    6 - 56                + |      3
    [2]             chr5    4 - 64                * |      4
    ------
    seqinfo(3 sequences): chr2 chr4 chr5

## Further information

- [Tutorial](https://biocpy.github.io/GenomicRanges/tutorial.html)
- [API documentation](https://biocpy.github.io/GenomicRanges/api/modules.html)
- [Bioc/GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
