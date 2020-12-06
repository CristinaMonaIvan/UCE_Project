


x <- read.table("indelsnv_seq_hg19_perpat.txt", sep="\t", quote="", head=T,check.names=F)

number <- x$mutations_per_mb
group <- x$type

median <- tapply(number, group, median)
group <- factor(group, levels=levels(group)[order(median, decreasing=T)])

y0 <- 10^-.5

number[number==0] <- y0

n10 <- number < 10
n100 <- 10 <= number & number < 100
n1000 <- 100 < number

color <- rep(NA, length(number))
color[n10] <- "darkseagreen1"
color[n100] <- "darkseagreen4"
color[n1000] <- "darkslategray"

m <- split(number, group)
col <- split(color, group)

tiff("fig3B.tif", width=6*300, height=6*300, res=300, compression="lzw")

par(font=2, font.axis=2, font.lab=2, mar=c(8,6,4,1)+.1)

plot(1, type="n", xlim=c(1, length(m)), ylim=c(y0, 1000),
		xlab="", ylab="",
		log="y", axes=F)

abline(h=c(y0, 1, 10, 100, 1000), col="grey90")

abline(v=seq(m), col="grey90")

for (i in seq(m)) {
	xc <- runif(length(m[[i]]), i-.3, i+.3)
	points(xc, m[[i]], col=col[[i]], pch=16,cex=1.5)
}

axis(1, at=seq(m), labels=names(m), lwd=2, las=2, cex.axis=1.5)

axis(2, at=c(y0, 1, 10, 100, 1000), labels=c(0, 1, 10, 100, 1000), lwd=2, las=1,cex.axis=1.5)

mtext("# UCE mutations per Mb per tumor", side=2, line=4, cex=1.6)

abline(h=c(10, 100), col=c("black", "darkslateblue"), lty="longdash")

box(lwd=2)

dev.off()
