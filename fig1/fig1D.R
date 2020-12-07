


source("source.R")


setwd("C:\\Users\\CAntonescu\\Desktop\\UCRpaper_data_analyses_writting\\ucr_coding_noncoding\\boxplots")



plegend <- function (p, prefix="P") {
	if (p<.0001) paste(prefix, "<0.0001", sep="") else paste(prefix, "=", formatC(p, digits=4, format="f"), sep="")
}


x <- read.table("forboxplot.txt", sep="\t", head=T, quote="\"", check.names=F, na="")

mir <- colnames(x)[4]

mirnames<-c(expression(bold(paste("# mutations per Mb per tumor"))))



names(mirnames) <- mir


status <- "Region"

r<-levels(x[,status])[c(1, 2,3,4)]

scale <- function (v) log(v+1)


for (m in mir) {
	l <- sapply(r, function (s) x[x[,status]%in%s, m], simplify=F)
	normal <- min(sapply(l, function (y) shapiro.test(y)$p.value)) > 0.05
	p <- (if (normal) t.test else wilcox.test)(l[[1]], l[[2]], paired=T)$p.value
	l <- lapply(l, scale)
	tiff(paste("UCRvsCoding_paired2",m,".tiff", sep="_"), width=11*300, height=2700, res=300, compression="lzw")
	par(oma=c(3,5,0,0),font=2)
	cex <- 2.5
	col=c("darkblue", "darkred", "darkmagenta", "palegreen4")
	r0 <- c(0, 1, 10,100, 1000)
	r1 <- scale(r0)
	boxplot(l, lwd=4, cex.axis=2.5, cex.main=2.5, cex=3, outline=T, ylim=range(r1), pch=c(18), main=plegend(p), axes=F, ylab="",
			boxcol=col,
			outcol=col,
			whiskcol=col,
			medcol=col,
			staplecol=col)
	for (i in seq(l[[1]])) lines(1:2, c(l[[1]][i], l[[2]][i]), lty=3)
	for (i in seq(l[[2]])) lines(2:3, c(l[[2]][i], l[[3]][i]), lty=3)
	for (i in seq(l[[3]])) lines(3:4, c(l[[3]][i], l[[4]][i]), lty=3)
	axis(2, at=r1,labels=r0, las=1, font.axis=2.5, cex.axis=cex, lwd=5)
	stripchart(l, vert=T, pch=c(17,18,19,20), cex=3, add=T, method="jitter", jitter=0, ylab="")
	mtext(mirnames[m], side=2, line=5, cex=2.5)
	mtext(c("coding", "UCE", "UCE\ncoding", "UCE\nnoncoding"), side=1, line=3, at=seq(l), cex=2.5)
	mtext("n=2449", side=1, line=6, at=2.47, cex=2.5)
	dev.off()
}


