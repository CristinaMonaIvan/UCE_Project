

x <- read.table("bed_type.txt", sep="\t", head=T, quote="\"", check.names=F, na="",stringsAsFactors=F)

m <- table(x$DISEASE, x$mut_type)

colnames(m) <- sub("_", ">", colnames(m))

col <- colors()[c(109,84,508,81,11,494,471,132)]

names(col) <- colnames(m)[order(colSums(m))]

m <- m[order(rownames(m), decreasing=T),order(colSums(m), decreasing=T)]

m <- apply(m, 1, function (v) 100*v/sum(v))


tiff("BarVariantClassPerCancerType.tif", width=9*300, height=4*300, res=300, compression="lzw")

par(font=2, font.axis=2, mar=c(5,7,7,3)+.1, mgp=c(9, 2, 0))

at <- seq(0, to=100, by=25)

b <- barplot(m, axes=F, xlim=range(at), cex.names=2.5, space=.4, col=col[rownames(m)], horiz=T, axisnames=F)

axis(1, cex.axis=2.7, lwd=6, at=at, labels=paste0(at, "%"))
mtext(colnames(m), side=2, line=2, cex=2.5, at=b, las=1)

legend(-23.5,5, rownames(m), col=col[rownames(m)], pch=15, bty="n", horiz=TRUE, x.intersp=.4, cex=2.1, xpd=T)

dev.off()

