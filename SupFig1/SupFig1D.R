

x <- read.table("mut_type.txt", sep="\t", head=T, quote="\"", check.names=F, na="",stringsAsFactors=F)

m <- table(x$type_pie, x$mut_type)

colnames(m) <- sub("_", ">", colnames(m))

col <- colors()[c(30, 109, 81, 508, 84, 11, 471, 494, 132)]

names(col) <- colnames(m)[order(colSums(m))]

m <- m[order(rownames(m), decreasing=T),order(colSums(m), decreasing=T)]

m <- apply(m, 1, function (v) 100*v/sum(v))


tiff("barSNV.tif", width=17*300, height=15*300, res=300, compression="lzw")

par(font=2, font.axis=2, mar=c(5,15,4,19)+.1, mgp=c(3, 2, 0))

at <- seq(0, to=100, by=25)

b <- barplot(m, axes=F, xlim=range(at), cex.names=2.7, space=.4, col=col[rownames(m)], horiz=T, axisnames=F)

axis(1, cex.axis=2.7, lwd=6, at=at, labels=paste0(at, "%"))
mtext(colnames(m), side=2, line=2, cex=2.7, at=b, las=1)

legend(105, 31, rownames(m), col=col[rownames(m)], pch=15, bty="n", cex=3, xpd=T)

dev.off()
