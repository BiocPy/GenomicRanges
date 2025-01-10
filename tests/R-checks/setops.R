x <- GRanges("chr1", IRanges(c(2, 9) , c(7, 19)), strand=c("+", "-"))
y <- GRanges("chr1", IRanges(5, 10), strand="-")

union(x, y)
union(x, y, ignore.strand=TRUE)

setdiff(x,y)
setdiff(x,y, ignore.strand=TRUE)

intersect(x,y)
intersect(x,y, ignore.strand=TRUE)
