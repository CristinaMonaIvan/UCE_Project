


source("sourcepie2D.R")

library(plotrix)

bracket <- function(x0, y0, x1, y0.top, y0.bottom, hook=5, lwd=2) {
	segments(x0, y0, x1, lwd=lwd)
	segments(x1, y0.top, y1=y0.bottom, lwd=lwd)
	segments(x1, y0.top, x1+hook, lwd=lwd)
	segments(x1, y0.bottom, x1+hook, lwd=lwd)
}

x <- read.table("number.txt", sep="\t",head=TRUE)

l <- x[,2]
names(l) <- paste(sub(" *$", "", x[,1]), " (n=", l, ")", sep="")

labels <- names(l)
labels[l<50] <- ""

n <- length(l)
cols <- rainbow(n)

r.adj <- rep(1, length(l))

#r.adj[10] <- 1.05

tiff("pie2dim.tif", width=14*300, height=12*300, res=300, compression="lzw")

par(font=2, mar=c(0,13.5,0,17)+.01)

pie(l, labels=labels, main="", col=cols, cex.main=3, cex=2, init.angle=68, clockwise=T)#, r.adj=r.adj)



zip <- function (v, n) {
	i <- c(t(suppressWarnings(matrix(seq(v), ncol=n))))
	v[i[seq(v)]]
}

ncol <- 4
nrow <- 2

legend(.1, -1.1, legend=zip(names(l)[l<50], nrow), pch=15, col=zip(cols[l<50], nrow),
		bty="n", ncol=ncol, cex=2, xpd=NA, x.intersp=.8, xjust=.5, text.width=.7)
dev.off()

