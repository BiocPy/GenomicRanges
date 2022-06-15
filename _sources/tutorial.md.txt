# Tutorial

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
            "seqname": ["chr1", "chr2", "chr3"],
            "start": [100, 115, 119],
            "end": [103, 116, 120],
        }
    )
)

hits = subject.nearest(query)
print(hits)
```

### Slice a GenomicRanges

You can slice a `GenomicRange` object using the subset (`[`) operator

```python
subject = GenomicRanges.fromUCSC(genome="hg38")

# first thousand features of the object
subset = subject[1:1000]
```
