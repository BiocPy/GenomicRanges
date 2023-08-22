[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)
[![PyPI-Server](https://img.shields.io/pypi/v/GenomicRanges.svg)](https://pypi.org/project/GenomicRanges/)
![Unit tests](https://github.com/BiocPy/GenomicRanges/actions/workflows/pypi-test.yml/badge.svg)

# GenomicRanges

Container class to represent genomic locations and support genomic analysis in Python similar to Bioconductor's [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).


## Install

Package is published to [PyPI](https://pypi.org/project/genomicranges/)

```shell
pip install genomicranges
```

## Usage

The package provide several ways to represent genomic annotations and intervals.

### Initialize a `GenomicRanges` object

#### From UCSC or GTF file

Methods are available to easily access UCSC genomes or load a genome annotation from GTF

```python
import genomicranges

gr = genomicranges.from_gtf(<PATH TO GTF>)
# OR
gr = genomicranges.from_ucsc(genome="hg19")
```
#### Pandas DataFrame

A common representation in Python is a pandas DataFrame for all tabular datasets. One can convert this into `GenomicRanges`. ***Intervals are inclusive on both ends.***

***Note: DataFrame must contain columns `seqnames`, `starts` and `ends` to represent genomic coordinates.***

```python
import genomicranges
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

gr = genomicranges.from_pandas(df)
```

### Interval Operations

Currently supports most commonly used [interval based operations](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).

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

Checkout the [documentation](https://biocpy.github.io/GenomicRanges/) for more usecases.


<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
