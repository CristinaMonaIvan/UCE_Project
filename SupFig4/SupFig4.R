

library(plotrix)

y <- read.table("table_all.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F)

w0 <- y[,c("change", "nucleotide")]
w0 <- w0[!duplicated(w0),]

w <- merge(data.frame(type=unique(y$type)), w0)

y <- merge(w, y, all.x=T)
y$count[is.na(y$count)] <- 0

col <- colors()[c(11,471,84,494,508,132)]
names(col) <- unique(y$change)


tiff("barplots.tiff", width=3*6*300, height=7*3*300, compression="lzw")

layout(matrix(1:24, ncol=3, byrow=T))

par(mar=c(23,6,4,0)+.1, cex=1, oma=c(0,3,3,0),
		font=2, font.axis=2, lwd=2)

for (i in unique(y$type)) {
	x <- y[y$type==i,]
	
	n <- x$count/sum(x$count)
	
	b <- barplot(n, col=col[x$change], axes=F, ylim=range(pretty(c(0,n))))
	
	axis(1, at=b, labels=F, lwd=5)
	axis(2, lwd=5, cex.axis=5, las=1, line=-3.5, at=pretty(c(0,n)))
	
	mtext(side=1, at=b, text=x$nucleotide, las=2, cex=3, line=c(2,8))
	
	y0 <- par("usr")[3]-3*diff(grconvertY(0:1, "inches", "user"))
	y1 <- par("usr")[3]-3.25*diff(grconvertY(0:1, "inches", "user"))
	x0 <- tapply(b, x$change, min)
	x1 <- tapply(b, x$change, max)
	
	rect(x0, y0, x1, y1, col=col, xpd=NA)
	
	mtext(side=1, text=names(col), line=20, at=(x0+x1)/2, cex=5)
	mtext(side=3, text=i, at=15, cex=7)
}

dev.off()