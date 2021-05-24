
library(RColorBrewer)

x <- read.table("forplot_2.txt", sep="\t", head=T, quote="", stringsAsFactors=T)

x <- x[order(x$NoList),]

colors <- c("Blues", "Greens", "YlOrBr", "Purples", "Reds", "Greys")

names(colors) <- levels(x$Data.base)

color <- function (x, lim, palette) palette[round(1+(length(palette)-1)*(x-lim[1])/diff(lim))]

pval.color <- function (p, pmin, pmax, pal.name) {
	pal <- colorRampPalette(brewer.pal(9, pal.name))(100)
	color(log(p), c(log(pmax), log(pmin)), pal)
}

pmin <- min(x$PValue)
pmax <- max(x$PValue)

col <- sapply(1:nrow(x), function (i) pval.color(x$PValue[i], pmin, pmax, colors[as.character(x$Data.base[i])]))


tiff("pathwaywithingenes.tif", width=28*300, height=18*300, res=300, compression="lzw")

layout(cbind(c(1,1,1), c(2,3,4), c(5,6,7)), widths=c(8.5,1))

par(cex.axis=1, mgp=c(3,1,0), mar=c(14,54,12,3)+.1, font=2, font.axis=2, oma=c(1,1,0,1), lwd=3, cex=1)

b <- barplot(
		x$NoList,
		names=paste0(x$Term, " (p=", format(x$PValue, digits=2), ")"),
		cex.names=3,
		las=1,
		horiz=T,
		space=0,
		main="",
		axes=F,
		lwd=5,
		cex.main=2.8,
		col=col,
		#xlim=range(pretty(x$No_genes_list))
		xlim=c(0,140)
)

axis(1, lwd=4, font=2, cex.axis=3, mgp=c(3,2,0), at=pretty(x$NoList))
mtext("Number of Genes", side=1, line=6, font=2, cex=3)

par(new=T)
par(font=2, font.axis=2)
v <- x$NoList/x$NoTotal
plot(v, b, type="l", axes=F, bty="n", xlab="", ylab="", xlim=range(pretty(v)))
axis(3, cex.axis=3, lwd=5, at=pretty(v), xpd=NA)
mtext("Ratio list/pathway", side=3, line=5, font=2, cex=3)

par(cex=1, mar=c(5.5, 5.5, 5.5, 5.5))

key <- seq(from=-log10(pmax), to=-log10(pmin), length=100)
for (c in names(colors))
{
	image(y=key,
			z=matrix(key, nrow=1),
			zlim=-log10(c(pmax, pmin)),
			col=colorRampPalette(brewer.pal(9, colors[c]))(100),
			useRaster=T, axes=F, xlab="", ylab="")
	axis(2, at=seq(from=2, to=8, by=1), las=1, cex.axis=2.8)
	mtext(expression(bold(paste("-log"["10"],"(p)"))), side=3, line=1, cex=2.8, xpd=T)
	mtext(c, side=1, line=1.5, cex=2.8, xpd=NA)
}

dev.off()
