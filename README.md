# GenomicRanges

Python equivalent to Bioconductor's [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) to represent genomic locations and support genomic analysis. It uses efficient structures already available in the Python/Pandas/numpy eco-system adds an familiar interfaces.


## Install

Package is deployed to [PyPI](https://pypi.org/project/genomicranges/)

```shell
pip install genomicranges
```

## Usage

The package provide several ways to represent genomic intervals

### Pandas DataFrame

A common representation in Python is a pandas DataFrame for all tabular datasets. One can convert this into `GenomicRanges`.

***Note: The DataFrame must contain columns `seqname`, `start` and `end` that represent chromosome and genomic coordinates.***

```python
from genomicranges import GenomicRanges

gr = GenomicRanges.fromPandas(<PANDAS DATA FRAME>)
```

### From UCSC or GTF file

Methods are available to easily access UCSC genomes or load a genome annotation from GTF

```python
from genomicranges import GenomicRanges

gr = GenomicRanges.fromGTF(<PATH TO GTF>)
# OR 
gr = GenomicRanges.fromUCSC(genome="hg19")
```

### Interval Operations

Currently supports [Nearest Genomic positions operation](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html#finding-the-nearest-genomic-position-in-granges-objects) in Bioconductor, but more coming soon.

```python
subject = GenomicRanges.fromUCSC(genome="hg38")

query = GenomicRanges.fromPandas(
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

For more use cases, checkout the [documentation](https://biocpy.github.io/GenomicRanges/)


<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
