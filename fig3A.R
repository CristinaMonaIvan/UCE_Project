



x <- read.table("indelsnv_seq_hg19_perpat.txt", sep="\t", quote="", head=T,check.names=F)

number <- sort(x$mutations_per_mb, decreasing=T)

y0 <- 10^-.5

number[number==0] <- y0

n10 <- number < 10
n100 <- 10 <= number & number < 100
n1000 <- 100 <= number

color <- rep(NA, length(number))
color[n10] <- "darkseagreen1"
color[n100] <- "darkseagreen4"
color[n1000] <- "darkslategray"

p <- function (n) 100*length(which(n))/length(number)

tiff("fig3A.tif", width=6*300, height=6*300, res=300, compression="lzw")

par(font=2, font.axis=2, font.lab=2, mar=c(8,6,4,1)+.1)

plot(seq(number), number, pch=16, col=color,
		xlab="number tumors", ylab="", cex.lab=1.7,
		log="y", axes=F)

abline(h=c(y0, 1, 10, 100, 1000), col="grey90")
abline(v=c(seq(from=0, to=2500, by=500)), col="grey90")

points(seq(number), number, pch=16, col=color,cex=1.5)

axis(1, lwd=2, cex.axis=1.5)
axis(2, at=c(y0, 1, 10, 100, 1000), labels=c(0, 1, 10, 100, 1000), lwd=2, las=1,cex.axis=1.5)

mtext("# UCE mutations per Mb per tumor", side=2, line=4, cex=1.6)

abline(h=c(10, 100), col=c("black", "darkslateblue"), lty="longdash")

text(length(number), c(10, 100, 1000), 
		paste0(formatC(c(p(n10), p(n100), p(n1000)), format="f", digits=1), "%"), 
		adj=c(.8, 1.5),cex=1.5)

box(lwd=2)

dev.off()
