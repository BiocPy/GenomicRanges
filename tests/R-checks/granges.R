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
GenomicRanges:::get_out_of_bound_index(x)
x_seqnames_id <- as.integer(seqnames(x))
x_seqlengths <- unname(seqlengths(x))
seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
seqlength_is_na <- is.na(x_seqlengths)
seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
which(seqlevel_has_bounds[x_seqnames_id] &
        (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))

idx <- GenomicRanges:::get_out_of_bound_index(x)

seqnames_id <- as.integer(seqnames(x))[idx]
new_end <- unname(seqlengths(x))[seqnames_id]
