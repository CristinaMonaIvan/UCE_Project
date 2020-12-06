

library(GenomicRanges)

x <- read.table("mut_all.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
y <- read.table("regions.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

q <- GRanges(x$CHROM, IRanges(x$POS, x$END))
s <- GRanges(y$chr, IRanges(y$start, y$end))

f <- findOverlaps(q, s, type="any")

m <- cbind(x[queryHits(f), c("CHROM", "POS")], y[subjectHits(f),"type",drop=F])
m <- m[!duplicated(m),]

write.table(m, file="mut_regions.txt", sep="\t", quote=F, row.names=F, na="")




