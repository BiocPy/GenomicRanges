# Tutorial

The package provide classes to represent genomic intervals and methods to perform interval based arithmetic operations.

**_Intervals are inclusive on both ends and starts at 1._**

The implementation details for the classes follow the Bioconductor's equivalent [R/GenomicRanges package](https://bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf).


# Construct a `GenomicRanges` object

To construct a **GenomicRanges** object, simply pass in the column representation as a dictionary. This dictionary must contain "seqnames", "starts", "ends" columns and optionally specify "strand". If "strand" column is not provided, "*" is used as the default value for each genomic interval.

### Pandas DataFrame

A common representation in Python is a pandas DataFrame for all tabular datasets. One can convert this DataFrame into `GenomicRanges` if it contains the necessary columns.

**_Note: The DataFrame must contain columns `seqnames`, `starts` and `ends` to represent genomic coordinates._**

```python
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

gr = GenomicRanges.from_pandas(df)
```

### From UCSC or GTF file

Methods are available to download, parse and access genomes from UCSC or load a genome annotation from GTF file.

```python
import genomicranges

gr = genomicranges.read_gtf(<PATH TO GTF>)
# OR
gr = genomicranges.read_ucsc(genome="hg19")
```

### from `IRanges` (Preferred way)

If you have all relevant information to create a ``GenomicRanges`` object

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
```

### Set sequence information

The package also provides a `SeqInfo` class to update or modify sequence information stored in the object. Read more about this class in [GenomeInfoDb package](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html).

```python
seq_obj = {
    "seqnames": ["chr1", "chr2", "chr3",],
    "seqlengths": range(100, 103),
    "is_circular": [random() < 0.5 for _ in range(3)],
    "genome": "hg19",
}

seq = SeqInfo(seq_obj)

gr.seq_info = seq
```

## Get and Set methods

Getters are available to access various properties.

```python
# access sequence names
gr.seqnames

# access all start positions
gr.start

# access annotation information if available
gr.seq_info

# compute and return the widths of each region
gr.width

# access metadata columns, everything other than genomic locations
gr.mcols
```

### Setters

All property based setters are `in_place` operations. Methods are available to get and set properties on GenomicRanges.

```python
gr.mcols = gr.mcols.set_column("score", range(1,6))

# or use an in-place operation
gr.mcols.set_column("score", range(1,6), in_place=True)
```

### Access any column

Aside from the default getters, `column` methods provides a way to quickly access any column in the object.

```python
gr.mcols.column("score")
```

### Access ranges

`ranges()` is a generic method to access only the genomic locations as dictionary, pandas `DataFrame` or something else. you can use any container representation based on a dictionary.

```python
gr.ranges
```

## Slice operations

You can slice a `GenomicRange` object using the subset (`[]`) operator. This operation accepts different slice input types, you can either specify a boolean vector, a `slice`` object, a list of indices, or row/column names to subset.

```python
# slice the first 3 rows
gr[:3]

# slice 1, 3 and 2nd rows
gr[[1,3,2]]
```

## Iterate over intervals

You can iterate over the intervals of a `GenomicRanges` object. `rowname` is None if the object does not have any row names.

```python
for rowname, row in gr:
    print(rowname, row)
```

## Intra-range transformations

For detailed description of these methods, refer to Bioconductor's [GenomicRanges documentation](https://bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf)

- **flank**: flank the intervals based on start or end or both.
- **shift**: shifts all the ranges specified by the shift argument.
- **resize**: resizes the ranges to the specified width where either the start, end, or center is used as an anchor
- **narrow**: narrows the ranges
- **promoters**: promoters generates promoter ranges for each range relative to the TSS.The promoter range is expanded around the TSS according to the upstream and downstream parameters.
- **restrict**: restricts the ranges to the interval(s) specified by the start and end arguments
- **trim**: trims out-of-bound ranges located on non-circular sequences whose length is not NA.

```python
# flank
flanked_gr = gr.flank(width=10, start=False, both=True)

# shift
shifted_gr = gr.shift(shift=10)

# resize
resized_gr = gr.resize(width=10, fix="end", ignore_strand=True)

# narrow
narrow_gr = gr.narrow(end=4, width=3)

# promoters
prom_gr = gr.promoters()

# restrict
restrict_gr = gr.restrict(start=114, end=140, keep_all_ranges=True)

# trim
trimmed_gr = gr.trim()
```

## Inter-range methods

- **range**: returns a new GenomicRanges object containing range bounds for each distinct (seqname, strand) pairing.
- **reduce**: returns a new GenomicRanges object containing reduced bounds for each distinct (seqname, strand) pairing.
- **gaps**: Finds gaps in the GenomicRanges object for each distinct (seqname, strand) pairing
- **disjoin**: Finds disjoint intervals across all locations for each distinct (seqname, strand) pairing.

```python
# range
range_gr = gr.range()

# reduce
reduced_gr = gr.reduce(min_gap_width=10, with_reverse_map=True)

# gaps
gapped_gr = gr.gaps(start=103) # OR
gapped_gr = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

# disjoin
disjoin_gr = gr.disjoin()
```

## Set operations on genomic ranges

- **union**: compute union of intervals across object
- **intersect**: compute intersection or finds overlapping intervals
- **setdiff**: compute set difference

```python
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

```python
# intersection
int_gr = g_src.intersect(g_tgt)

# set diff
diff_gr = g_src.setdiff(g_tgt)

# union
union_gr = g_src.union(g_tgt)
```

## Compute over bins

### Summary stats for column

one can use Pandas for this

```python
pd.Series(gr.mcols.get_column("score")).describe()
```

### `binned_average`

Compute binned average for different positions

```python
bins = pd.DataFrame({"seqnames": ["chr1"], "starts": [101], "ends": [109],})

bins_gr = GenomicRanges.from_pandas(bins)

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
binned_avg_gr
```

now you might wonder how can I generate these ***bins***?

### Generate tiles or bins from `GenomicRanges`

- **tile**: Splits each genomic region by n (number of regions) or by width (maximum width of each tile)
- **sliding_windows**: Generates sliding windows within each range, by width and step.

```python
# tiles
tiles = gr.tile(n=2)

# slidingwindows
tiles = gr.sliding_windows(width=10)
```

### Generate tiles from Genome

`tile_genome` returns a set of genomic regions that form a partitioning of the specified genome.

```python
seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

tiles = GenomicRanges.tile_genome(seqlengths=seqlengths, n=10)
```

### Coverage

Computes number of ranges that overlap for each position in the range.

```python
res_vector = gr.coverage(shift=10, width=5)
```

## Overlap based methods

- **find_overlaps**: find overlaps between two `GenomicRanges` object
- **count_overlaps**: count overlaps between two `GenomicRanges` object
- **subset_by_overlaps**: subset a `GenomicRanges` object if it overlaps with the ranges in the query

```python
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

# findOverlaps
res = subject.find_overlaps(query, query_type="within")

# countOverlaps
res = subject.count_overlaps(query)

# subsetByOverlaps
res = subject.subset_by_overlaps(query)
```

## Search operations

- **nearest**: Performs nearest neighbor search along any direction (both upstream and downstream)
- **follow**: Performs nearest neighbor search only along downstream
- **precede**: Performs nearest neighbor search only along upstream

```python
find_regions = GenomicRanges(
    seqnames= ["chr1", "chr2", "chr3"],
    ranges=IRanges([200, 105, 1190],[203, 106, 1200]),
)

query_hits = gr.nearest(find_regions)

query_hits = gr.precede(find_regions)

query_hits = gr.follow(find_regions)
```

## Comparison, rank and order operations

- **match**: Element wise comparison to find exact match intervals.
- **order**: Get the order of indices for sorting.
- **sort**: Sort the GenomicRanges object.
- **rank**: for each interval identifies its position is a sorted order

```python
# match
query_hits = gr.match(gr[2:5])

# order
order = gr.order()

# sort
sorted_gr = gr.sort()

# rank
rank = gr.rank()
```

## Combine `GenomicRanges` objects by rows

Use the `combine` generic from [biocutils](https://github.com/BiocPy/generics) to concatenate multiple GenomicRanges objects.

```python
from biocutils.combine import combine
combined_gr = combine(gr, gr1, gr2, ...)
```

or use the `combine`` method,

```python
combined_gr = gr.combine(gr1, gr2, gr3, ...)
```

## Misc operations

- **invert_strand**: flip the strand for each interval
- **sample**: randomly choose ***k*** intervals

```python
# invert strand
inv_gr = gr.invert_strand()

# sample
samp_gr = gr.sample(k=4)
```

# Construct a `GenomicRangesList` object.

Just as it sounds, a `GenomicRangesList` is a named-list like object.

If you are wondering why you need this class, a `GenomicRanges` object lets us specify multiple
genomic elements, usually where the genes start and end. Genes are themselves made of many sub
regions, e.g. exons. `GenomicRangesList` allows us to represent this nested structure.

Currently, this class is limited in functionality, purely a read-only class with basic accessors.

***Note: This is a work in progress and the functionality is limited.***

```python
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
grl
```


## Properties

```python
grl.start
grl.width
```

## Combine `GenomicRangeslist` object

Similar to the combine function from GenomicRanges,

```python

grla = GenomicRangesList(ranges=[a], names=["a"])
grlb = GenomicRangesList(ranges=[b, a], names=["b", "c"])

# or use the combine generic
from biocutils.combine import combine
cgrl = combine(grla, grlb)
```

and that's all for now! Check back later for more updates.
