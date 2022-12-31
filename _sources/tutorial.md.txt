# Tutorial

The package provide classes to represent genomic locations and methods to perform interval based operations. ***Intervals are inclusive on both ends.***

For detailed description of these methods, checkout [GenomicRanges documentation](https://bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf)

# Construct a `GenomicRanges` object
## Pandas DataFrame

A common representation in Python is a pandas DataFrame for all tabular datasets. One can convert this DataFrame into `GenomicRanges`.

***Note: The DataFrame must contain columns `seqnames`, `starts` and `ends` to represent genomic coordinates.***

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

gr = genomicranges.fromPandas(df)
```

## From UCSC or GTF file

Methods are available to access UCSC genomes or load a genome annotation from GTF

```python
import genomicranges

gr = genomicranges.readGTF(<PATH TO GTF>)
# OR 
gr = genomicranges.readUCSC(genome="hg19")
```

## Set sequence information

```python
seq_obj = {
    "seqnames": ["chr1", "chr2", "chr3",],
    "seqlengths": range(100, 103),
    "isCircular": [random() < 0.5 for _ in range(3)],
    "genome": "hg19",
}

seq = SeqInfo(seq_obj)

gr.seqInfo = seq
```

# Get and Set methods

Getters are available to access various properties.

```python
# access sequence names
gr.seqnames

# access all start positions
gr.starts

# access annotation information if available
gr.seqInfo

# compute and return the widths of each region
gr.width

# score if available
gr.score

# access metadata columns, everything other than genomic locations
gr.mcols()
```

## Setters

Set properties in class

```python
gr.score = [<NEW ARRAY OF VALUES>]
```

## Access any column

Aside from the default getters, `column` methods provides a way to quicky access any column in the object.

```python
gr.column("seqnames")

gr.column("score")
```

## Access ranges 

`ranges()` is a generic method to access only the genomic locations as dictionary, pandas `DataFrame` or something else. you can use any container representation based on a dictionary.

```python
# default to dict
gr.ranges()

# as pandas DataFrame
gr.ranges(returnType=pd.DataFrame)
```

`granges()` method returns a new `GenomicRanges` object of just the genomic locations

```python
gr.granges()
```

# Slice operations

You can slice a `GenomicRange` object using the subset (`[`) operator.

```python
# slice the first 5 rows
gr[:5, :]

# slice 1,3 and 5th rows
gr[[1,3,5], :]
```

# Iterate over rows

```python
for index, row in gr:
    print(index, row)
```


# Intra-range transformations

For detailed description of these methods, checkout [GenomicRanges documentation](https://bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf)

- flank: flank the intervals based on start or end or both. 
- shift: shifts all the ranges specified by the shift argument.
- resize: resizes the ranges to the specified width where either the start, end, or center is used as an anchor
- narrow: narrows the ranges
- promoters: promoters generates promoter ranges for each range relative to the TSS.The promoter range is expanded around the TSS according to the upstream and downstream parameters.
- restrict: restricts the ranges to the interval(s) specified by the start and end arguments
- trim: trims out-of-bound ranges located on non-circular sequences whose length is not NA.


```python
# flank
flanked_gr = gr.flank(width=10, start=False, both=True)

# shift
shifted_gr = gr.shift(shift=10)

# resize
resized_gr = gr.resize(width=10, fix="end", ignoreStrand=True)

# narrow
narrow_gr = gr.narrow(end=4, width=3)

# promoters
prom_gr = gr.promoters()

# restrict
restrict_gr = gr.restrict(start=114, end=140, keepAllRanges=True)

# trim
trimmed_gr = gr.trim()
```

# Inter-range methods

- range: returns a new GenomicRanges object containing range bounds for each distinct (seqname, strand) pairing.
- reduce: returns a new GenomicRanges object containing reduced bounds for each distinct (seqname, strand) pairing.
- gaps: Finds gaps in the GenomicRanges object for each distinct (seqname, strand) pairing
- disjoin: Finds disjoint intervals across all locations for each distinct (seqname, strand) pairing.
- isDisjoint: Is the object contain disjoint intervals for each distinct (seqname, strand) pairing?

```python
# range
range_gr = gr.range()

# reduce
reduced_gr = gr.reduce(minGapwidth=10, withRevMap=True)

# gaps
gapped_gr = gr.gaps(start=103) # OR 
gapped_gr = gr.gaps(end={"chr1": 120, "chr2": 120, "chr3": 120})

# disjoin
disjoin_gr = gr.disjoin()

# isDisjoint
isdisjoin = gr.isDisjoint()
```

# Set operations on genomic ranges

- union: compute union of intervals across object
- intersect: compute intersection or finds overlapping intervals
- setdiff: compute set difference 

```python
df_src = pd.DataFrame(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3", "chr2"],
        "starts": [101, 102, 103, 104, 109],
        "ends": [112, 103, 128, 134, 111],
        "strand": ["*", "-", "*", "+", "-"],
        "score": range(0, 5),
        "GC": [random() for _ in range(5)],
    }
)

g_src = genomicranges.fromPandas(df_src)

df_tgt = pd.DataFrame(
    {
        "seqnames": ["chr1","chr2","chr2","chr2","chr1","chr1","chr3","chr3","chr3","chr3"],
        "starts": range(101, 111),
        "ends": range(121, 131),
        "strand": ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

g_tgt = genomicranges.fromPandas(df_tgt)
```

```python
# union
union_gr = g_src.union(g_tgt)

int_gr = g_src.intersect(g_tgt)

diff_gr = g_src.setdiff(g_tgt)
```

# Compute over bins

## Summary stats for column

one can use Pandas for this

```python
pd.Series(gr.column("score")).describe()
```

## `binnedAverage`

compute binned average for different positions

```python
bins = pd.DataFrame({"seqnames": ["chr1"], "starts": [101], "ends": [109],})

bins_gr = genomicranges.fromPandas(bins)

subject = pd.DataFrame(
    {
        "seqnames": ["chr1","chr2","chr2","chr2","chr1","chr1","chr3","chr3","chr3","chr3"],
        "starts": range(101, 111),
        "ends": range(121, 131),
        "strand": ["*", "-", "-", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

subject_gr = genomicranges.fromPandas(subject)


# Compute binned average
binned_avg_gr = g_tgt.binnedAverage(bins=bins_gr, scorename="score", outname="binned_score")
```

now you might wonder how can I generate these bin?

## Generate tiles or bins from `GenomicRanges`

- `tile`: splits each genomic region by n (number of regions) or by width (maximum width of each tile)
- `slidingWindows`: Generates sliding windows within each range, by width and step.

```python
# tiles
tiles = gr.tile(n=2)

# slidingwindows
tiles = gr.slidingWindows(width=10)
```

## Generate tiles from Genome

`tileGenome` returns a set of genomic regions that form a partitioning of the specified genome.


```python
seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

tiles = genomicranges.tileGenome(seqlengths=seqlengths, n=10)
```

## Coverage

computes number of ranges that overlap for each position in the range.

```python
res_vector = gr.coverage(shift=10, width=5)
```

# Overlap based methods

- findOverlaps: find overlaps between two GenomicRanges object
- countOverlaps: count overlaps between two GenomicRanges object
- subsetByOverlaps: subset a GenomicRanges object if it overlaps with the ranges in the query

```python
df_subject = pd.DataFrame(
    {
        "seqnames": [
            "chr1",
            "chr2",
            "chr2",
            "chr2",
            "chr1",
            "chr1",
            "chr3",
            "chr3",
            "chr3",
            "chr3",
        ],
        "starts": range(1, 11),
        "ends": [10] * 10,
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"],
        "score": range(0, 10),
        "GC": [random() for _ in range(10)],
    }
)

subject = genomicranges.fromPandas(df_subject)

df_query = pd.DataFrame(
    {"seqnames": ["chr2",], "starts": [4], "ends": [6], "strand": ["+"]}
)

query = genomicranges.fromPandas(df_query)
```

```python
# findOverlaps
res = subject.findOverlaps(query, queryType="within")

# countOverlaps
res = subject.countOverlaps(query)

# subsetByOverlaps
res = subject.subsetByOverlaps(query)
```

# Search operations

- nearest: Performs nearest neighbor search along any direction (both upstream and downstream)
- follow: Performs nearest neighbor search only along downstream
- precede: Performs nearest neighbor search only along upstream
- distanceToNearest: calculate distance to nearest location


```python

find_regions = genomicranges.fromPandas(
    pd.DataFrame(
        {
            "seqnames": ["chr1", "chr2", "chr3"],
            "starts": [200, 105, 1190],
            "ends": [203, 106, 1200],
        }
    )
)

query_hits = gr.nearest(find_regions)

query_hits = gr.precede(test_gr)

query_hits = gr.follow(test_gr)

query_hits = gr.distanceToNearest(test_gr)
```

# Comparison, rank and order operations

- duplicated: if any of the ranges are duplicated
- match: Element wise comparison to find exact match intervals.
- isUnsorted: if the object is not sorted
- order: Get the order of indices for sorting.
- sort: Sort the GenomicRanges object.
- rank: for each interval identifies its position is a sorted order

```python
# duplicated
query_hits = gr.duplicated()

# match
query_hits = gr.match(gr[2:5, :])

# is unsorted?
result = gr.isUnsorted()

# order
order = gr.order()

# sort
sorted_gr = gr.sort()

# rank
rank = gr.rank()
```

# Concatenate GenomicRanges objects by rows

concatenate any number of GenomicRanges objects

```python
concat_gr = gr.concat([gr1, gr2, gr3, ...])
```

# Misc operations

- invertStrand: flip the strand for each interval
- sample: randomly choose k intervals

```python
# invert strand
inv_gr = gr.invertStrand()

# sample
samp_gr = gr.sample(k=4)
```