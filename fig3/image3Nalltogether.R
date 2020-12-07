

source("source_splitscale.R")


y <- read.table("table_all.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F)

x <- aggregate(y[,c("count"),drop=F], y[,c("change", "nucleotide")], sum)
x <- x[order(x$change, x$nucleotide),]

ticks <- list(seq(from=0, to=.03, by=.005), c(.05, 0.055), c(.08, .085))
ranges <- t(sapply(ticks, range))
scale <- split.scale(ranges, heights=sapply(ticks, length)-1, space=1)

n <- scale(x$count/sum(x$count))

col <- colors()[c(11,471,84,494,508,132)]
names(col) <- unique(x$change)

tiff("HighHyper.tiff", width=12*300, height=6*300, compression="lzw")

par(mar=c(32,9,4,0)+.1,
		font=2, font.axis=2, lwd=2)

b <- barplot(n, col=col[x$change], axes=F, ylim=range(scale(ranges)))

axis(1, at=b, labels=F, lwd=4)

for (r in ticks) axis(2, lwd=4, cex.axis=5, las=1, line=-5, at=scale(r), labels=format(r))

mtext(side=2, text="Relative Contribution", cex=5, line=10)
mtext(side=1, at=b, text=x$nucleotide, las=2, cex=5, line=c(2,12))

y0 <- par("usr")[3]-4.45*diff(grconvertY(0:1, "inches", "user"))
y1 <- par("usr")[3]-4.7*diff(grconvertY(0:1, "inches", "user"))
x0 <- tapply(b, x$change, min)
x1 <- tapply(b, x$change, max)

rect(x0, y0, x1, y1, col=col, xpd=NA)

mtext(side=1, text=names(col), line=27, at=(x0+x1)/2, cex=5)
mtext(side=1, text="96 trinucleotide changes", line=33, cex=5)

dev.off()

