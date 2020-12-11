



x <- read.table("mut_regions.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)


y <- x[,c(1,2,3)]
y <- y[duplicated(y),]
y <- y[!duplicated(y),]
z <- merge(x, y)

write.table(y, file="morethanoneVAR.txt", sep="\t", quote=F, row.names=F, na="")
