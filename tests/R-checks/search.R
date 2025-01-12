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


gr = GRanges(
  seqnames=c("chr1","chr2","chr3","chr2","chr3"),
  ranges=IRanges(c(101, 102, 103, 104, 105), width=c(11, 21, 25, 30, 5)),
  strand=c("*", "-", "*", "+", "-")
)


x = IRanges(c(1, 5, -2, 0, 14), width=c(10, 5, 6, 12, 4))

seqlengths <- c(chr1=60, chr2=20, chr3=25)

## Create 5 tiles:
tiles <- tileGenome(seqlengths, ntile=5)
tiles

x <- GRanges(c(A="chr1:1-50", B="chr1:40-110", C="chrX:1-500"))
y <- GRanges(c("chr1:21-25", "chr1:38-150"))
z <- subtract(x, y)
z

ignore.strand <- FALSE
minoverlap <- 1
y_red <- reduce(y, ignore.strand=ignore.strand)
hits <- findOverlaps(x, y_red, minoverlap=minoverlap,
                            ignore.strand=ignore.strand)
hits_obj <- extractList(y, as(hits, "IntegerList"))
setdiff(x[1], hits_obj[[1]])
ans <- psetdiff(x, extractList(y, as(hits, "IntegerList")))
mcols(ans) <- mcols(x)
setNames(ans, names(x))
