

plegend <- function (p, prefix="P") {
	if (p<.0001) paste(prefix, "<0.0001", sep="") else paste(prefix, "=", formatC(p, digits=4, format="f"), sep="")
}


y <- read.table("count.txt", sep="\t", head=T, quote="\"", check.names=F, na="")

number <- c(y$count_CDS,
		y$count_IGR,
		y$count_Intron,
		y$count_ncRNA,
		y$count_UTR)

class <- rep(c("CDS", "IGR", "Intron", "ncRNA", "3p&5p"), each=nrow(y))

x <- data.frame(number=number, class=class)

mir <- colnames(x)[1]
mirnames<-c("# UCE mutations")
names(mirnames) <- mir

x[,2]<-factor(x[,2])


status <- "class"

r<-levels(x[,status])

scale <- function (v) log(v+1)

d2<-NULL

for (m in mir) {
	l <- sapply(r, function (s) x[x[,status]%in%strsplit(s, "\\.")[[1]], m], simplify=F)
	if (any(sapply(l, function (y) length(which(y!=0))<3))) next()
	normal <- min(sapply(l, function (y) shapiro.test(y)$p.value)) > 0.05
	p <- if (normal) aov.pvalue(aov.test(l)) else kruskal.test(l)$p.value
	l <- lapply(l, scale)
	tiff(paste(m,"VariantClassification.tiff", sep="_"), width=15*300, height=2700, res=300, compression="lzw")
	par(oma=c(3,5,0,0),font=2)
	cex <- 2.8
	#r1 <- pretty(do.call(c, l))
	r0 <- c(0, 1, 10,100, 1000)
	r1 <- scale(r0)
	b <- boxplot(l, lwd=4, cex.axis=2.8, cex.main=2.8, outline=T, ylim=range(r1), main=plegend(p), col=c("darkred", "magenta", "slategrey", "darkgreen", "darkblue"), axes=F, ylab="")
	axis(2, at=r1, labels=r0, las=1, font.axis=2, cex.axis=cex, lwd=5)
	mtext(mirnames[m], side=2, line=6, cex=3)
	mtext(c("3p&5p","CDS", "IGR", "Intron", "ncRNA"), side=1, line=3.5, at=seq(l), cex=3)
	dev.off()
}


