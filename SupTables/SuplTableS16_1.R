




x <- read.table("intervals_r.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

x <- x[!is.na(x$start)&!is.na(x$end),]

x <- x[order(as.numeric(sub("X", 23, sub("Y", 24, x$chr))), x$start, -x$end),]

write.table(x, file="intervals_sorted.txt", sep="\t", quote=F, row.names=F, na="")
