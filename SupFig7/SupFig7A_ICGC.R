

library(gplots)

x <- read.table("table.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
z <- read.table("table_type.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

x <- x[order(as.numeric(sub("X", "23", sub("chr", "", x$chr))), x$start),]

chr <- tapply(1:nrow(x), x$chr, range)
chr <- chr[unique(x$chr)]
chr.mid <- sapply(chr, mean)
chr.sep <- sapply(chr, "[", 1)[-1]

y <- t(as.matrix(x[,5:26]))

tiff("allUCRevents.tiff", width=2000, height=4000, res=300, compression="lzw")

layout(matrix(1:2), heights=c(12,1))

par(mar=c(7,4,4,4)+.1)

n <- 5

col <- c("white", colorpanel(n, "green", "purple"))

image(1:nrow(y), 1:ncol(y), pmin(y,n), col=col, zlim=c(0,n), axes=F, xlab="", ylab="")

axis(1, at=1:nrow(y), labels=paste0(rownames(y), "(", z$Patients, ")"), las=2)

odd <- seq(from=1, to=length(chr), by=2)
even <- seq(from=2, to=length(chr), by=2)
axis(2, at=chr.mid[odd], labels=names(chr.mid)[odd], las=1, lwd.ticks=0, mgp=c(3,.5,0))
axis(4, at=chr.mid[even], labels=names(chr.mid)[even], las=1, lwd.ticks=0, mgp=c(3,.5,0))

abline(v=2:(ncol(y)-1)-.5, lwd=1,lty=2)
abline(h=chr.sep, lwd=1,lty=2)

par(mar=c(3,4,1,4))
image(x=0:n, y=1, z=matrix(0:n), col=col, axes=F, xlab="", ylab="")
axis(1, at=0:n)
box()

dev.off()
