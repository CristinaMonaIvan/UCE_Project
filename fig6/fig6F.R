

library(reshape2)
library(RColorBrewer)

x <- read.table("table.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)


p <- melt(x[,c(1, grep("Pvalue", colnames(x)))], id.vars=1, variable.name="Gene", value.name="Pvalue")
p$Gene <- sub(" Pvalue", "", p$Gene)

fch <- melt(x[,c(1, grep("FCH", colnames(x)))], id.vars=1, variable.name="Gene", value.name="FCH")
fch$Gene <- sub(" FCH", "", fch$Gene)

y <- merge(p, fch)

blue <- colorRampPalette(brewer.pal(9, "Blues"))(100)[20:90]
red <- colorRampPalette(brewer.pal(9, "Reds"))(100)[20:90]

color <- function (z, lim, palette) palette[round(1+(length(palette)-1)*(z-lim[1])/diff(lim))]

plim <- range(-log10(y$Pvalue))
pcolor <- function (p, palette) color(-log10(p), plim, palette)
col <- function (p, fch) ifelse(fch>0, pcolor(p, red), pcolor(p, blue))

i <- seq(unique(y$Gene))
names(i) <- unique(y$Gene)

j <- c(7,6,5,4,3,2,1)
names(j) <- unique(y$Study)


lab <- paste(x$Study, x$Samples, sep="\n")
names(lab) <- x$Study

tiff("fig6F.tiff", width=300*9, height=300*7, res=300, compression="lzw")

layout(matrix(c(1, 1, 1, 2, 3, 4), ncol=2), widths=c(8, 1))

par(mar=c(5, 23.3, 4, 3)+.1, font=2, font.lab=2, font.axis=2, lwd=2, cex=1.027)

plot(i[y$Gene], j[y$Study], xlim=c(min(i)-.2, max(i)), pch=16, cex=1.75*abs(y$FCH), col=col(y$Pvalue, y$FCH),
		axes=F, xlab="", ylab="", xpd=NA)

axis(1, at=i, labels=names(i), lwd=2)
axis(2, at=j, labels=lab[names(j)], las=2, cex=2, lwd=2)

par(mar=c(1.67,2,3,2.5))

plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=F)
legend(.5, .5, legend=c(1, 2, 3), pch=16, pt.cex=1.5*c(1,2,3),
		xjust=.5, yjust=.5, cex=1.2, bty="n", xpd=NA, y.intersp=2, adj=c(-1, .5))
text(.45, 1.07, "Fold Change", xpd=NA)

par(mar=c(2.5,2,2,2.5))

key <- seq(from=min(plim), to=max(plim), length=100)
image(y=key, z=matrix(key, nrow=1), zlim=range(plim), col=red,  useRaster=T, axes=F, xlab="", ylab="")
mtext(expression(bold(paste("-log"["10"],"(p)"))), side=3, line=.5, cex=1.2, xpd=T)
axis(2, las=1, cex.axis=1, font.axis=2)

par(mar=c(3,2,1.5,2.5))

key <- seq(from=min(plim), to=max(plim), length=100)
image(y=key, z=matrix(key, nrow=1), zlim=range(plim), col=blue,  useRaster=T, axes=F, xlab="", ylab="")
mtext(expression(bold(paste("-log"["10"],"(p)"))), side=3, line=.5, cex=1.2, xpd=T)
axis(2, las=1, cex.axis=1, font.axis=2)

dev.off()
