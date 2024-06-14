---
file_format: mystnb
kernelspec:
  name: python
---

# `GenomicRanges`: Genomic analysis

`GenomicRanges` is a Python package designed to handle genomic locations and facilitate genomic analysis. It is similar to Bioconductor's [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) and uses the [IRanges](https://github.com/BiocPy/IRanges) package under the hood to manage and provide interval-based arithmetic operations.

An `IRanges` holds a **start** position and a **width**, and is typically used to represent coordinates along a genomic sequence. The interpretation of the **start** position depends on the application; for sequences, the **start** is usually a 1-based position, but other use cases may allow zero or even negative values, e.g., circular genomes. `IRanges` uses [nested containment lists](https://github.com/pyranges/ncls) under the hood to perform fast overlap and search-based operations.

The package provides a `GenomicRanges` class to specify multiple genomic elements, typically where genes start and end. Genes are themselves made of many subregions, such as exons, and a `GenomicRangesList` enables the representation of this nested structure.

Moreover, the package also provides a `SeqInfo` class to update or modify sequence information stored in the object. Learn more about this in the [GenomeInfoDb package](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html).

The `GenomicRanges` class is designed to seamlessly operate with upstream packages like `RangeSummarizedExperiment` or `SingleCellExperiment` representations, providing consistent and stable functionality.

:::{note}
These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.
:::

## Installation

To get started, install the package from [PyPI](https://pypi.org/project/genomicranges/)

```bash
pip install genomicranges
```

# Construct a `GenomicRanges` object

We support multiple ways to initialize a `GenomicRanges` object.

## From Bioinformatic file formats

### From `biobear`

Although the parsing capabilities in this package are limited, the [biobear](https://github.com/wheretrue/biobear) library is designed for reading and searching various bioinformatics file formats, including FASTA, FASTQ, VCF, BAM, and GFF, or from an object store like S3. Users can esily convert these representations to `GenomicRanges` (or [read more here](https://www.wheretrue.dev/docs/exon/biobear/genomicranges-integration)):

```python
from genomicranges import GenomicRanges
import biobear as bb

session = bb.new_session()

df = session.read_gtf_file("path/to/test.gtf").to_polars()
df = df.rename({"seqname": "seqnames", "start": "starts", "end": "ends"})

gg = GenomicRanges.from_polars(df)

# do stuff w/ a genomic ranges
print(len(gg), len(df))
```

### From UCSC or GTF file

You can also import genomes from UCSC or load a genome annotation from a GTF file. This requires installation of additional packages **pandas** and **joblib** to parse and extract various attributes from the gtf file.

```python
import genomicranges

# gr = genomicranges.read_gtf(<PATH TO GTF>)

# OR

human_gr = genomicranges.read_ucsc(genome="hg19")
print(human_gr)
```


## Preferred way

To construct a `GenomicRanges` object, we need to provide sequence information and genomic coordinates. This is achieved through the combination of the `seqnames` and `ranges` parameters. Additionally, you have the option to specify the `strand`, represented as a list of "+" (or 1) for the forward strand, "-" (or -1) for the reverse strand, or "*" (or 0) if the strand is unknown. You can also provide a NumPy vector that utilizes either the string or numeric representation to specify the `strand`. Optionally, you can use the `mcols` parameter to provide additional metadata about each genomic region.

```{code-cell}
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

:::{note}
The input for `mcols` is expected to be a `BiocFrame` object and will be converted to a `BiocFrame` in case a pandas `DataFrame` is supplied.
:::

## Pandas `DataFrame`

If your genomic coordinates are represented as a pandas `DataFrame`, convert this into `GenomicRanges` if it contains the necessary columns.

::: {important}
The `DataFrame` must contain columns `seqnames`, `starts` and `ends` to represent genomic coordinates. The rest of the columns are considered metadata and will be available in the `mcols` slot of the `GenomicRanges` object.
:::

```{code-cell}
from genomicranges import GenomicRanges
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

gr_from_df = GenomicRanges.from_pandas(df)
print(gr_from_df)
```

## Polars `DataFrame`

Similarly, To initialize from a polars `DataFrame`:

```{code-cell}
from genomicranges import GenomicRanges
import polars as pl
from random import random

df = pl.DataFrame(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3", "chr2"],
        "starts": [101, 102, 103, 104, 109],
        "ends": [112, 103, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

gr = GenomicRanges.from_polars(df)
print(gr)
```

# Sequence information

The package also provides a `SeqInfo` class to update or modify sequence information stored in the object. Learn more about this in the [GenomeInfoDb package](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html).

```{code-cell}
from genomicranges import SeqInfo

seq = SeqInfo(
    seqnames = ["chr1", "chr2", "chr3"],
    seqlengths = [110, 112, 118],
    is_circular = [True, True, False],
    genome = "hg19",
)
gr_with_seq = gr.set_seqinfo(seq)
print(gr_with_seq)
```

## Getters/Setters

Getters are available to access various attributes using either the property notation or functional style.

```{code-cell}
# access sequence names
print("seqnames (as property): ", gr.seqnames)
print("seqnames (functional style): ", gr.get_seqnames())

# access all start positions
print("start positions: ", gr.start)

# access annotation information if available
gr.seqinfo

# compute and return the widths of each region
print("width of each region: ", gr.get_width())
# or gr.width

# access mcols
print(gr.mcols)
```

### Setters

:::{important}
All property-based setters are `in_place` operations, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.
:::

```{code-cell}
modified_mcols = gr.mcols.set_column("score", range(1,6))
modified_gr = gr.set_mcols(modified_mcols)
print(modified_gr)
```

or use an in-place operation:

```{code-cell}
gr.mcols.set_column("score", range(1,6), in_place=True)
print(gr.mcols)
```

### Access ranges

`get_ranges()` is a generic method to access only the genomic coordinates:

```{code-cell}
# or gr.get_ranges()
print(gr.ranges)
```

# Subset operations

You can subset a `GenomicRange` object using the subset (`[]`) operator. This operation accepts different slice input types, such as a boolean vector, a `slice` object, a list of indices, or names (if available) to subset.

```{code-cell}
# get the first 3 regions
gr[:3]

# get 1, 3 and 2nd rows
# note: the order is retained in the result
print(gr[[1,3,2]])
```

## Iterate over ranges

You can iterate over the regions of a `GenomicRanges` object. `name` is `None` if the object does not contain any names. To iterate over the first two ranges:

```{code-cell}
for name, row in gr[:2]:
    print(name, row)
```

# Intra-range transformations

For detailed description of these methods, refer to either the Bioconductor's or BiocPy's documentation.

- **flank**: Flank the intervals based on **start** or **end** or **both**.
- **shift**: Shifts all the ranges specified by the **shift** argument.
- **resize**: Resizes the ranges to the specified width where either the **start**, **end**, or **center** is used as an anchor.
- **narrow**: Narrows the ranges.
- **promoters**: Promoters generates promoter ranges for each range relative to the TSS. The promoter range is expanded around the TSS according to the **upstream** and **downstream** parameters.
- **restrict**: Restricts the ranges to the interval(s) specified by the **start** and **end** arguments.
- **trim**: Trims out-of-bound ranges located on non-circular sequences whose length is not `NA`.

```{code-cell}
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

# flank
flanked_gr = gr.flank(width=10, start=False, both=True)

# shift
shifted_gr = gr.shift(shift=10)

# resize
resized_gr = gr.resize(width=10, fix="end", ignore_strand=True)

# narrow
narrow_gr = gr.narrow(end=1, width=1)

# promoters
prom_gr = gr.promoters()

# restrict
restrict_gr = gr.restrict(start=114, end=140, keep_all_ranges=True)

# trim
trimmed_gr = gr.trim()

print("GenomicRanges after the trim operation:")
print(trimmed_gr)
```

# Inter-range methods

- **range**: Returns a new `GenomicRanges` object containing range bounds for each distinct (seqname, strand) pair.
- **reduce**: returns a new `GenomicRanges` object containing reduced bounds for each distinct (seqname, strand) pair.
- **gaps**: Finds gaps in the `GenomicRanges` object for each distinct (seqname, strand) pair.
- **disjoin**: Finds disjoint intervals across all locations for each distinct (seqname, strand) pair.

```{code-cell}
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

# range
range_gr = gr.range()

# reduce
reduced_gr = gr.reduce(min_gap_width=3, with_reverse_map=True)

# gaps
gapped_gr = gr.gaps(start=103)  # OR
gapped_gr = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

# disjoin
disjoin_gr = gr.disjoin()

print("GenomicRanges with the disjoint ranges:")
print(disjoin_gr)
```

# Set operations

- **union**: Compute the `union` of intervals across objects.
- **intersect**: Compute the `intersection` or finds overlapping intervals.
- **setdiff**: Compute `set difference`.

```{code-cell}
#| code-fold: true
#| code-summary: "Show the code"

g_src = GenomicRanges(
    seqnames = ["chr1", "chr2", "chr1", "chr3", "chr2"],
    ranges = IRanges(start =[101, 102, 103, 104, 109], width=[112, 103, 128, 134, 111]),
    strand = ["*", "-", "*", "+", "-"]
)

g_tgt = GenomicRanges(
    seqnames = ["chr1","chr2","chr2","chr2","chr1","chr1","chr3","chr3","chr3","chr3"],
    ranges = IRanges(start =range(101, 111), width=range(121, 131)),
    strand = ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"]
)
```

```{code-cell}
# intersection
int_gr = g_src.intersect(g_tgt)

# set diff
diff_gr = g_src.setdiff(g_tgt)

# union
union_gr = g_src.union(g_tgt)

print("GenomicRanges after the union operation:")
print(union_gr)
```

# Compute over bins

## Summary stats for column

Use `Pandas` to compute summary statistics for a column:

```{code-cell}
pd.Series(gr.mcols.get_column("score")).describe()
```

With a bit more magic, render a histogram using **matplotlib**:

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

_ = plt.hist(gr.mcols.get_column("score"), bins="auto")
plt.title("'score' histogram with 'auto' bins")
plt.show()
```

Not the prettiest plot but it works.

## Binned average

Compute binned average for a set of query **bins**:

```{code-cell}
from iranges import IRanges
bins_gr = GenomicRanges(seqnames=["chr1"], ranges=IRanges([101], [109]))

subject = GenomicRanges(
    seqnames= ["chr1","chr2","chr2","chr2","chr1","chr1","chr3","chr3","chr3","chr3"],
    ranges=IRanges(range(101, 111), range(121, 131)),
    strand= ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame({
        "score": range(0, 10),
    })
)

# Compute binned average
binned_avg_gr = subject.binned_average(bins=bins_gr, scorename="score", outname="binned_score")
print(binned_avg_gr)
```

::: {tip}
Now you might wonder how can I generate these ***bins***?
:::

# Generate tiles or bins

- **tile**: Splits each genomic region by **n** (number of regions) or by **width** (maximum width of each tile).
- **sliding_windows**: Generates sliding windows within each range, by **width** and **step**.

```{code-cell}
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

# tiles
tiles = gr.tile(n=2)

# slidingwindows
tiles = gr.sliding_windows(width=10)
print(tiles)
```

## Generate tiles from genome

`tile_genome` returns a set of genomic regions that form a partitioning of the specified genome.

```{code-cell}
seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

tiles = GenomicRanges.tile_genome(seqlengths=seqlengths, n=10)
print(tiles)
```

# Coverage

Computes number of ranges that overlap for each position.

```{code-cell}
import rich

res_vector = gr.coverage()
rich.print(res_vector)
```

Lets see what the coverage looks like, now with **seaborn**:

```{code-cell}
import seaborn as sns
vector = res_vector["chr1"]
sns.lineplot(data=pd.DataFrame({
    "position": [i for i in range(len(vector))],
    "coverage":vector
}), x ="position", y="coverage")
```

I guess that looks ok. :) but someone can help make this visualization better. (something that ports `plotRanges` from R)

# Overlap based methods

- **find_overlaps**: Find overlaps between two `GenomicRanges` objects.
- **count_overlaps**: Count overlaps between two `GenomicRanges` objects.
- **subset_by_overlaps**: Subset a `GenomicRanges` object if it overlaps with the ranges in the query.

```{code-cell}
subject = GenomicRanges(
    seqnames= ["chr1","chr2","chr2","chr2","chr1","chr1","chr3","chr3","chr3","chr3"],
    ranges=IRanges(range(101, 111), range(121, 131)),
    strand= ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
    mcols=BiocFrame({
        "score": range(0, 10),
    })
)

df_query = pd.DataFrame(
    {"seqnames": ["chr2",], "starts": [4], "ends": [6], "strand": ["+"]}
)

query = GenomicRanges.from_pandas(df_query)

# find Overlaps
res = subject.find_overlaps(query, query_type="within")

# count Overlaps
res = subject.count_overlaps(query)

# subset by Overlaps
res = subject.subset_by_overlaps(query)

print(res)
```

# Search operations

- **nearest**: Performs nearest neighbor search along any direction (both upstream and downstream).
- **follow**: Performs nearest neighbor search only along downstream.
- **precede**: Performs nearest neighbor search only along upstream.

```{code-cell}
find_regions = GenomicRanges(
    seqnames= ["chr1", "chr2", "chr3"],
    ranges=IRanges([200, 105, 1190],[203, 106, 1200]),
)

query_hits = gr.nearest(find_regions)

query_hits = gr.precede(find_regions)

query_hits = gr.follow(find_regions)

print(query_hits)
```

::: {note}
Similar to `IRanges` operations, these methods typically return a list of indices from `subject` for each interval in `query`.
:::

# Comparison, rank and order operations

- **match**: Element-wise comparison to find exact match intervals.
- **order**: Get the order of indices for sorting.
- **sort**: Sort the `GenomicRanges` object.
- **rank**: For each interval identifies its position is a sorted order.

```{code-cell}
# match
query_hits = gr.match(gr[2:5])
print("matches: ", query_hits)

# order
order = gr.order()
print("order:", order)

# sort
sorted_gr = gr.sort()
print("sorted:", sorted_gr)

# rank
rank = gr.rank()
print("rank:", rank)
```

# Combine `GenomicRanges` objects by rows

Use the `combine` generic from [biocutils](https://github.com/BiocPy/generics) to concatenate multiple `GenomicRanges` objects.

```{code-cell}

from biocutils.combine import combine
a = GenomicRanges(
    seqnames=["chr1", "chr2", "chr1", "chr3"],
    ranges=IRanges([1, 3, 2, 4], [10, 30, 50, 60]),
    strand=["-", "+", "*", "+"],
    mcols=BiocFrame({"score": [1, 2, 3, 4]}),
)

b = GenomicRanges(
    seqnames=["chr2", "chr4", "chr5"],
    ranges=IRanges([3, 6, 4], [30, 50, 60]),
    strand=["-", "+", "*"],
    mcols=BiocFrame({"score": [2, 3, 4]}),
)

combined = combine(a,b)
print(combined)
```

# Misc operations

- **invert_strand**: flip the strand for each interval
- **sample**: randomly choose ***k*** intervals

```{code-cell}
# invert strand
inv_gr = gr.invert_strand()

# sample
samp_gr = gr.sample(k=4)
```

# `GenomicRangesList` class

Just as it sounds, a `GenomicRangesList` is a named-list like object.

If you are wondering why you need this class, a `GenomicRanges` object enables the
specification of multiple genomic elements, usually where genes start and end.
Genes, in turn, consist of various subregions, such as exons.
The `GenomicRangesList` allows us to represent this nested structure.

As of now, this class has limited functionality, serving as a read-only class with basic accessors.

```{code-cell}

from genomicranges import GenomicRangesList
a = GenomicRanges(
    seqnames=["chr1", "chr2", "chr1", "chr3"],
    ranges=IRanges([1, 3, 2, 4], [10, 30, 50, 60]),
    strand=["-", "+", "*", "+"],
    mcols=BiocFrame({"score": [1, 2, 3, 4]}),
)

b = GenomicRanges(
    seqnames=["chr2", "chr4", "chr5"],
    ranges=IRanges([3, 6, 4], [30, 50, 60]),
    strand=["-", "+", "*"],
    mcols=BiocFrame({"score": [2, 3, 4]}),
)

grl = GenomicRangesList(ranges=[a,b], names=["gene1", "gene2"])
print(grl)
```


## Properties

```{code-cell}
grl.start
grl.width
```

## Combine `GenomicRangeslist` object

Similar to the combine function from `GenomicRanges`,

```{code-cell}
grla = GenomicRangesList(ranges=[a], names=["a"])
grlb = GenomicRangesList(ranges=[b, a], names=["b", "c"])

# or use the combine generic
from biocutils.combine import combine
cgrl = combine(grla, grlb)
```

The functionality in `GenomicRangesLlist` is limited to read-only and a few methods. Updates are expected to be made as more features become available.

## Empty ranges

Both of these classes can also contain no range information, and they tend to be useful when incorporates into larger data structures but do not contain any data themselves.

To create an empty `GenomicRanges` object:

```{code-cell}
empty_gr = GenomicRanges.empty()

print(empty_gr)
```

Similarly, an empty `GenomicRangesList` can be created:

```{code-cell}
empty_grl = GenomicRangesList.empty(n=100)

print(empty_grl)
```

----

## Futher reading

- Check out [the reference documentation](https://biocpy.github.io/GenomicRanges/) for more details.
- Visit Bioconductor's [**GenomicRanges**](https://bioconductor.org/packages/GenomicRanges) package.
