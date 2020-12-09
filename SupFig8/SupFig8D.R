


plegend <- function (p, prefix="P") {
	if (p<.0001) paste(prefix, "<0.0001", sep="") else paste(prefix, "=", formatC(p, digits=4, format="f"), sep="")
}


x <- read.table("CLLclinicalinformation_NoPreviousTreatm.txt", sep="\t", head=T, quote="\"", check.names=F, na="")

mir <- colnames(x)[34]
mirnames<-c(expression(bold("# UCE mutations per Mb per tumor")))
names(mirnames) <- mir

title <- ""

status <- "percentageZAP70_grp"

r<-levels(x[,status])

scale <- function (v) log(v+1)

d2<-NULL

for (m in mir) {
	l <- sapply(r, function (s) x[x[,status]%in%strsplit(s, "\\.")[[1]], m], simplify=F)
	if (any(sapply(l, function (y) length(which(y!=0))<3))) next()
	normal <- min(sapply(l, function (y) shapiro.test(y)$p.value)) > 0.05
	p <- (if (normal) t.test else wilcox.test)(l[[1]], l[[2]])$p.value
	l <- lapply(l, scale)
	tiff(paste(m,"ZAP70_2.tiff", sep="_"), width=9*300, height=2700, res=300, compression="lzw")
	par(oma=c(3,5,0,0),font=2)
	cex <- 2.8
	r0 <- c(10, 20, 40)
	r1 <- scale(r0)
	b <- boxplot(l, lwd=4, cex.axis=2.8, cex.main=2.8, outline=T, ylim=range(r1), main=plegend(p), col=c("lightblue", "firebrick1"), axes=F, ylab="")
	axis(2, at=r1, labels=r0, las=1, font.axis=2.8, cex.axis=cex, lwd=5)
	mtext(mirnames[m], side=2, line=6, cex=2.8)
	mtext(c("ZAP70-", "ZAP70+"), side=1, line=3, at=seq(l), cex=2.8)
	mtext(paste("n=", sapply(l, length), sep=""), side=1, line=6, at=seq(l), cex=2.8)
	dev.off()
}




