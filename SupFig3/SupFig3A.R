
x <- read.table("indelsnv_seq_hg19_perpat.txt", sep="\t", quote="\"", head=T, check.names=F)

number <- x$mutations_per_mb
group <- x$type

group<-factor(group)

m<-mean(number)
s<-sd(number)
max<-max(number)
p <- c(m, m+s, m+2*s)

col<-rainbow(length(levels(group)))

# by median
median <- tapply(number, group, median)
o <- order(-median, levels(group), decreasing=T)
group <- factor(group, levels=levels(group)[o])
col <- col[o]
labels <- levels(group)

scale <- function (v) log(1+v)

tiff("boxplot.tiff", width=14*300, height=10*300, res=300, compression="lzw")

par(mar=c(3,13,1,2)+.1, font=2)

boxplot(scale(number)~group, 
		horizontal=T,
		cex=1.5,
		axes=F,
		outpch=16,
		boxfill=col,
		outcol=col
)

at <- c(0, 10, 100, 1000)

axis(1, at=scale(at), labels=at, font.axis=2, lwd=2, cex.axis=1.7)
axis(2, at=seq(labels), labels=rep("", length(labels)), las=1, font.axis=2, lwd=2,cex.axis=1.5)

mtext(paste0(labels," (n=", table(group), ")"), side=2, line=1, at=seq(labels), las=1, cex=1.5)

mtext("# UCR mutations per Mb per tumor", side=1, line=3, las=1, cex=1.7)

abline(v=scale(p), lty=c("dotted", "dotdash","longdash"), lwd=2)

legend(4.8, 4, 
		lty=c("dotted", "dotdash","longdash", "blank"),
		lwd=2,
		legend=paste0(c("mean=", "mean+SD=", "mean+2SD=", "max="), 
				formatC(c(p,max), format="f", digits=0)),
		cex=1.5,
		bty="n",
		xpd=NA)

dev.off()

