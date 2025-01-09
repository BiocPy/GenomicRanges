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

ir <- IRanges(1:10, width=10:1)
restrict(ir, start=4, end=8, keep.all.ranges=TRUE)
