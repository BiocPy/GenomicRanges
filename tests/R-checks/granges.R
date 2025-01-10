library(GenomicRanges)
gr <- GRanges(
  seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
  ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
  strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score=1:10,
  GC=seq(1, 0, length=10)
)
gr

flank(gr, width=10)
flank(gr, width=10, start=FALSE)
flank(gr, width=10, both=TRUE)
flank(gr, width=10, start=FALSE, both=TRUE)

resize(gr, width=10)
resize(gr, width=10, fix="end")
resize(gr, width=10, ignore.strand=TRUE)
resize(gr, width=10, fix="end", ignore.strand=TRUE)
resize(gr, width=10, fix="center", ignore.strand=TRUE)
resize(gr, width=10, fix="center")
resize(gr, width=11)
resize(gr, width=11, fix="center")

narrow(gr, start=2, end=3)
narrow(gr, start=2)
narrow(gr, start=2, width=3)
narrow(gr, end=2)
narrow(gr, width=3, end=4)


shift(gr, shift=10)
promoters(gr)
paste(as.vector(start(promoters(gr))), collapse=",")
paste(as.vector(end(promoters(gr))), collapse=",")

terminators(gr)
paste(as.vector(start(terminators(gr))), collapse=",")
paste(as.vector(end(terminators(gr))), collapse=",")

restrict_gr = restrict(gr, start=3, end=7)
paste(as.vector(start(restrict_gr)), collapse=",")
paste(as.vector(end(restrict_gr)), collapse=",")

restrict_gr = restrict(gr, start=4, end=8, keep.all.ranges=TRUE)
restrict_gr
paste(as.vector(start(restrict_gr)), collapse=",")
paste(as.vector(end(restrict_gr)), collapse=",")

restrict_gr = restrict(gr, start=1200)
paste(as.vector(start(restrict_gr)), collapse=",")
paste(as.vector(end(restrict_gr)), collapse=",")

st <- structure(c(4,5), names = c("chr1", "chr2"))
en <-  structure(c(8,9), names = c("chr2", "chr3"))

restrict_gr = restrict(gr, start=st, end=en)
paste(as.vector(start(restrict_gr)), collapse=",")
paste(as.vector(end(restrict_gr)), collapse=",")

gr = GRanges(
  seqnames=c(
    "chr1",
    "chr2",
    "chr3",
    "chr2",
    "chr3"
  ),
  ranges=IRanges(101:105, width=c(11, 21, 25, 30, 5)),
  strand=c("*", "-", "*", "+", "-")
)

seq_obj = Seqinfo(
  seqnames=c("chr1", "chr2", "chr3"),
  seqlengths=c(110, 112, 118),
  isCircular=c(TRUE, TRUE, FALSE),
  genome="hg19"
)
seq_obj
x <- gr
seqinfo(x) <- seq_obj
trim(x)

## reduce
reduce(gr)
reduce(gr, ignore.strand=TRUE)
reduce(gr, min.gapwidth=10)
reduce(gr, min.gapwidth=10, with.revmap=TRUE)

gr2 <- GRanges(
  seqnames=c(
    "chr1_gl123",
    "chr2",
    "chr3",
    "chr2",
    "chr3"),
  ranges=IRanges(101:105, width=c(11, 21, 25, 30, 5)),
  strand=c("*", "-", "*", "+", "-")
)
reduce(gr2)
reduce(gr2, ignore.strand=TRUE)


range(gr)
range(gr, ignore.strand=TRUE)

gaps(gr)
gaps(gr, ignore.strand=TRUE)
gaps(gr, start=5)
gaps(gr, start=103)
gaps(gr, end=c(chr1= 120, chr2= 120, chr3= 120))
gaps(gr, start=5, end=10)

disjoin(gr)
disjoin(gr, ignore.strand=TRUE)

isDisjoint(gr)
isDisjoint(gr, ignore.strand=TRUE)

disjointBins(gr)
disjointBins(gr, ignore.strand=TRUE)

coverage(gr)
coverage(gr, shift = 10)
coverage(gr, shift = 10, width=5)
