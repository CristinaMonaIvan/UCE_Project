x <- read.table("severalF006_columns.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

f <- grepl("F004|F005", x$FILTER)

n <- aggregate(list(f=f), x[,1:2], any)
n <- n[n$f, -3]

y <- merge(x, n)
y <- y[order(y$CHROM, y$POS),]

write.table(y, file="severalF006_2.txt", sep="\t", quote=F, row.names=F, na="")
