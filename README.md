[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)
[![PyPI-Server](https://img.shields.io/pypi/v/GenomicRanges.svg)](https://pypi.org/project/GenomicRanges/)
![Unit tests](https://github.com/BiocPy/GenomicRanges/actions/workflows/pypi-test.yml/badge.svg)

# GenomicRanges

GenomicRanges provides container classes designed to represent genomic locations and support genomic analysis. It is similar to Bioconductor's [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).

**_Intervals are inclusive on both ends and starts at 1._**

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

### Pandas DataFrame

A common representation in Python is a pandas `DataFrame` for all tabular datasets. `DataFrame` must contain columns "seqnames", "starts", and "ends" to represent genomic intervals. Here's an example:

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
print(gr)
```

    ## output
    GenomicRanges with 5 intervals & 2 metadata columns                                
    ┏━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┓
    ┃ row_names ┃ seqnames <list> ┃ starts <list> ┃ ends <list> ┃ strand <list> ┃ score <list> ┃ GC <list>           ┃
    ┡━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━┩
    │ 0         │ chr1            │ 101           │ 112         │ *             │ 0            │ 0.22617584001235103 │
    │ 1         │ chr2            │ 102           │ 103         │ -             │ 1            │ 0.25464256182466394 │
    │ ...       │ ...             │ ...           │ ...         │ ...           │ ...          │ ...                 │
    │ 4         │ chr2            │ 109           │ 111         │ -             │ 4            │ 0.5414168889911801  │
    └───────────┴─────────────────┴───────────────┴─────────────┴───────────────┴──────────────┴─────────────────────┘

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

## `GenomicRangesList`

Just as it sounds, a `GenomicRangesList` is a named-list like object. If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple genomic elements, usually where the genes start and end. Genes are themselves made of many sub-regions, e.g. exons. `GenomicRangesList` allows us to represent this nested structure.

**Currently, this class is limited in functionality.**

To construct a GenomicRangesList

```python
gr1 = GenomicRanges(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3"],
        "starts": [1, 3, 2, 4],
        "ends": [10, 30, 50, 60],
        "strand": ["-", "+", "*", "+"],
        "score": [1, 2, 3, 4],
    }
)

gr2 = GenomicRanges(
    {
        "seqnames": ["chr2", "chr4", "chr5"],
        "starts": [3, 6, 4],
        "ends": [30, 50, 60],
        "strand": ["-", "+", "*"],
        "score": [2, 3, 4],
    }
)

grl = GenomicRangesList(ranges=[gr1, gr2], names=["gene1", "gene2"])
print(grl)
```

    ## output
    GenomicRangesList with 2 genomic elements                                                                          
                                                                                                                    
    Name: gene1                                                                                            
                GenomicRanges with 4 intervals & 1 metadata columns                                           
    ┏━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┓                                  
    ┃ seqnames <list> ┃ starts <list> ┃ ends <list> ┃ strand <list> ┃ score <list> ┃      
    ┡━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━┩                                  
    │ chr1            │ 1             │ 10          │ -             │ 1            │                                  
    │ chr2            │ 3             │ 30          │ +             │ 2            │                                  
    │ chr3            │ 4             │ 60          │ +             │ 4            │                                  
    └─────────────────┴───────────────┴─────────────┴───────────────┴──────────────┘                                  
                                                                                                                    
    Name: gene2                                                                                            
                GenomicRanges with 3 intervals & 1 metadata columns                                           
    ┏━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┓                                  
    ┃ seqnames <list> ┃ starts <list> ┃ ends <list> ┃ strand <list> ┃ score <list> ┃      
    ┡━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━┩                                  
    │ chr2            │ 3             │ 30          │ -             │ 2            │                                  
    │ chr4            │ 6             │ 50          │ +             │ 3            │                                  
    │ chr5            │ 4             │ 60          │ *             │ 4            │                                  
    └─────────────────┴───────────────┴─────────────┴───────────────┴──────────────┘             

## Further information

- [Tutorial](https://biocpy.github.io/GenomicRanges/tutorial.html)
- [API documentation](https://biocpy.github.io/GenomicRanges/api/modules.html)
- [Bioc/GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.1.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
