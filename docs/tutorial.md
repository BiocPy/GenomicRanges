# Tutorial

The package provide classes to represent genomic locations and methods to perform interval based operations. **_Intervals are inclusive on both ends and starts at 1._**

For detailed description of these methods, checkout [R/GenomicRanges documentation](https://bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf).


# Construct a `GenomicRanges` object

To construct a **GenomicRanges** object, simply pass in the column representation as a dictionary. This dictionary must contain "seqnames", "starts", "ends" columns and optionally specify "strand". If "strand" column is not provided, "*" is used as the default value for each genomic interval.

## Pandas DataFrame

A common representation in Python is a pandas DataFrame for all tabular datasets. One can convert this DataFrame into `GenomicRanges`.

**_Note: The DataFrame must contain columns `seqnames`, `starts` and `ends` to represent genomic coordinates._**

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

## From UCSC or GTF file

Methods are available to access UCSC genomes or load a genome annotation from GTF

```python
import genomicranges

gr = genomicranges.read_gtf(<PATH TO GTF>)
# OR
gr = genomicranges.read_ucsc(genome="hg19")
```

## From a dictionary

```python
gr = GenomicRanges(
    {
        "seqnames": ["chr1", "chr2", "chr3"],
        "starts": [100, 115, 119],
        "ends": [103, 116, 120],
    }
)
```


## Set sequence information

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

# Get and Set methods

Getters are available to access various properties.

```python
# access sequence names
gr.seqnames

# access all start positions
gr.starts

# access annotation information if available
gr.seq_info

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

Aside from the default getters, `column` methods provides a way to quickly access any column in the object.

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
gr.ranges(return_type=pd.DataFrame)
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

# Inter-range methods

- **range**: returns a new GenomicRanges object containing range bounds for each distinct (seqname, strand) pairing.
- **reduce**: returns a new GenomicRanges object containing reduced bounds for each distinct (seqname, strand) pairing.
- **gaps**: Finds gaps in the GenomicRanges object for each distinct (seqname, strand) pairing
- **disjoin**: Finds disjoint intervals across all locations for each distinct (seqname, strand) pairing.
- **is_disjoint**: Is the object contain disjoint intervals for each distinct (seqname, strand) pairing?

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

# is Disjoint?
isdisjoin = gr.is_disjoint()
```

# Set operations on genomic ranges

- **union**: compute union of intervals across object
- **intersect**: compute intersection or finds overlapping intervals
- **setdiff**: compute set difference

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

g_src = genomicranges.from_pandas(df_src)

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

g_tgt = genomicranges.from_pandas(df_tgt)
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

## `binned_average`

Compute binned average for different positions

```python
bins = pd.DataFrame({"seqnames": ["chr1"], "starts": [101], "ends": [109],})

bins_gr = genomicranges.from_pandas(bins)

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

subject_gr = genomicranges.from_pandas(subject)


# Compute binned average
binned_avg_gr = g_tgt.binned_average(bins=bins_gr, scorename="score", outname="binned_score")
```

now you might wonder how can I generate these ***bins***?

## Generate tiles or bins from `GenomicRanges`

- **tile**: Splits each genomic region by n (number of regions) or by width (maximum width of each tile)
- **sliding_windows**: Generates sliding windows within each range, by width and step.

```python
# tiles
tiles = gr.tile(n=2)

# slidingwindows
tiles = gr.sliding_windows(width=10)
```

## Generate tiles from Genome

`tile_genome` returns a set of genomic regions that form a partitioning of the specified genome.

```python
seqlengths = {"chr1": 100, "chr2": 75, "chr3": 200}

tiles = genomicranges.tile_genome(seqlengths=seqlengths, n=10)
```

## Coverage

Computes number of ranges that overlap for each position in the range.

```python
res_vector = gr.coverage(shift=10, width=5)
```

# Overlap based methods

- **find_overlaps**: find overlaps between two `GenomicRanges` object
- **count_overlaps**: count overlaps between two `GenomicRanges` object
- **subset_by_overlaps**: subset a `GenomicRanges` object if it overlaps with the ranges in the query

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

subject = genomicranges.from_pandas(df_subject)

df_query = pd.DataFrame(
    {"seqnames": ["chr2",], "starts": [4], "ends": [6], "strand": ["+"]}
)

query = genomicranges.from_pandas(df_query)
```

```python
# findOverlaps
res = subject.find_overlaps(query, queryType="within")

# countOverlaps
res = subject.count_overlaps(query)

# subsetByOverlaps
res = subject.subset_by_overlaps(query)
```

# Search operations

- **nearest**: Performs nearest neighbor search along any direction (both upstream and downstream)
- **follow**: Performs nearest neighbor search only along downstream
- **precede**: Performs nearest neighbor search only along upstream
- **distance_to_nearest**: calculate distance to nearest location

```python

find_regions = genomicranges.from_pandas(
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

query_hits = gr.distance_to_nearest(test_gr)
```

# Comparison, rank and order operations

- **duplicated**: if any of the ranges are duplicated
- **match**: Element wise comparison to find exact match intervals.
- **is_unsorted**: if the object is not sorted
- **order**: Get the order of indices for sorting.
- **sort**: Sort the GenomicRanges object.
- **rank**: for each interval identifies its position is a sorted order

```python
# duplicated
query_hits = gr.duplicated()

# match
query_hits = gr.match(gr[2:5, :])

# is unsorted?
result = gr.is_unsorted()

# order
order = gr.order()

# sort
sorted_gr = gr.sort()

# rank
rank = gr.rank()
```

# Concatenate `GenomicRanges` objects by rows

Concatenate any number of GenomicRanges objects,

```python
concat_gr = gr.concat([gr1, gr2, gr3, ...])
```

# Misc operations

- **invert_strand**: flip the strand for each interval
- **sample**: randomly choose k intervals

```python
# invert strand
inv_gr = gr.invert_strand()

# sample
samp_gr = gr.sample(k=4)
```

# Construct a `GenomicRangesList` object.

***Note: This is a work in progress and the functionality is very limited.***

```python
a = GenomicRanges(
    {
        "seqnames": ["chr1", "chr2", "chr1", "chr3"],
        "starts": [1, 3, 2, 4],
        "ends": [10, 30, 50, 60],
        "strand": ["-", "+", "*", "+"],
        "score": [1, 2, 3, 4],
    }
)

b = GenomicRanges(
    {
        "seqnames": ["chr2", "chr4", "chr5"],
        "starts": [3, 6, 4],
        "ends": [30, 50, 60],
        "strand": ["-", "+", "*"],
        "score": [2, 3, 4],
    }
)

grl = GenomicRangesList(gene1=a, gene2=b)

# access normal properties same as GenomicRanges

grl.start
grl.width

```
