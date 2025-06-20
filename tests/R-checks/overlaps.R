library(GenomicRanges)
gr <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges=IRanges(1:10, width=10:1, names=head(letters,10)),
  strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score=1:10,
  GC=seq(1, 0, length=10)
)
gr

## GRangesList object:
gr1 <- GRanges(seqnames="chr2", ranges=IRanges(4:3, 6),
               strand="+", score=5:4, GC=0.45)
gr1

findOverlaps(gr1, gr)
findOverlaps(gr1, gr, ignore.strand=TRUE)
findOverlaps(gr, gr1, ignore.strand=TRUE)

findOverlaps(gr1, gr, type = "within")

x = GRanges("chr1", IRanges(c(2, 9), width=c(7, 19)), strand=c("+", "-"))
y = GRanges("chr1", IRanges(5, width=10), strand=c("*"))

findOverlaps(x, y)
findOverlaps(y, x)

countOverlaps(gr1, gr)
