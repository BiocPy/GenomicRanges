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

nearest(gr1, gr)
nearest(gr, gr1)
nearest(gr1, gr, select="all")
nearest(gr1, gr, ignore.strand=TRUE)
nearest(gr1, gr, select="all", ignore.strand=TRUE)

precede(gr1, gr)
precede(gr, gr1)
precede(gr1, gr, select="all")

follow(gr1, gr)
follow(gr, gr1)
follow(gr1, gr, select="all")
