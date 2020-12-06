

library(plotrix)

y <- read.table("table_all.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F)

w0 <- y[,c("change", "nucleotide")]
w0 <- w0[!duplicated(w0),]

w <- merge(data.frame(type=unique(y$type)), w0)

y <- merge(w, y, all.x=T)
y$count[is.na(y$count)] <- 0

z <- aggregate(count~type, y, sum)
colnames(z)[2] <- "total"

y <- merge(y, z)
y$ratio <- y$count/y$total

col <- rainbow(22)
names(col) <- unique(y$type)

x <- y[,c("change", "nucleotide")]
x <- x[!duplicated(x),]
x <- x[order(x$change, x$nucleotide),]
x$coord <- 1:nrow(x)

y <- merge(y, x)

spread <- .15

tiff("pointplot.tiff", width=20*300, height=7*300, res=300, compression="lzw")

par(font=2, font.axis=2, lwd=3, xaxs="i", yaxs="i",
		mar=c(5, 4, 2, 2)+.1)

plot(y$coord+runif(length(y$coord),-spread, spread),
		y$ratio,
		xlim=range(c(.5, x$coord+.5)),
		ylim=range(c(-.005, pretty(y$ratio))),
		pch=19,
		col=col[y$type],
		cex=1.2,
		xlab="",
		ylab="",
		bty="n",
		axes=F)

axis(1, at=x$coord, labels=F)
axis(2, las=1, line=0)

mtext(side=1, at=x$coord, text=x$nucleotide, las=2, cex=1, line=.7)

cl <- tapply(x$coord, x$change, min)
cr <- tapply(x$coord, x$change, max)
cm <- (cl+cr)/2

mtext(side=1, at=cm, text=names(cm), cex=1.5, line=3.7)

rect(cl, -.035, cr, -.033, col="black", xpd=NA)

abline(v=seq(from=min(x$coord)-.5, to=max(x$coord)+.5, by=1), lwd=1, lty=2)
abline(v=cl[-1]-.5, lwd=2)

box()


dev.off()

tiff("legend.tif", width=1.5*300, height=5*300, res=300, compression="lzw")

par(mar=c(0,1,0,1))

plot.new()

legend(.5, .5, xjust=.5, yjust=.5, legend=names(col), col=col, pch=19, xpd=NA, bty="n", text.font=2)
#legend(.4, .4, xjust=.5, yjust=.5, pch=15, col=col, legend=c("A", "C", "G", "T"), pt.cex=5.5, text.font=2, cex=2, bty="n", xpd=NA)

dev.off()
