


source("mozaicplot.R")

plegend <- function (p, prefix="P") {
	if (p<.001) paste(prefix, "<0.001", sep="") else paste(prefix, "=", formatC(p, digits=3, format="f"), sep="")
}

z <- read.table("fisher50species.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F, na="", check.names=F)


tiff("species.tiff", width=2100, height=2700, res=300, compression="lzw")


par(font=2, mar=c(5,7.5,4,2)+.1, cex.main=2.2)

t <- t(matrix(as.numeric((z[1,c(2,1,4,3)])), nrow=2))

rownames(t) <- paste0(c("Spp<50", "Spp=50"), "\n(n=", rowSums(t), ")")
colnames(t) <- c("", "")

rect.label <- t(apply(t, 1, function (v) formatC(v/sum(v), digits=2, format="f")))

mosaicplot(t, off=c(1,0), xlab="", ylab="", main="", , cex=2.2, color=c("firebrick1", "lightblue"),
		rect.label=rect.label, rect.cex=2.2, rect.font=2)

legend(-.38, .6, legend=c("UCR MUT", "no UCR\n MUT"), pch=15, col=c("firebrick1", "lightblue"), bty="n", xpd=T, cex=2.2)


legend(-.45, .95, legend=plegend(z$pvalue), xpd=T, bty="n", cex=2.2)

dev.off()


