

plegend <- function (p)(paste0("=",formatC(p, format="e", digits=1)))

x <- read.table("indelsnv_seq_hg19_allcoding_perpat.txt", sep="\t", quote="\"", head=T, check.names=F)
p <- read.table("pvalues.txt", sep="\t", quote="\"", head=T, check.names=F, stringsAsFactors=F)

number <- c(x$mutations_per_mb_Coding,
         x$mutations_per_mb_UCR)


group <- factor(rep(x$type, 2))
class <- rep(c("Coding", "UCR"), each=nrow(x))

stat <- NULL
for (i in unique(class)) {
	m<-mean(number[class==i])
	s<-sd(number[class==i])
	max<-max(number[class==i])
	stat <- rbind(stat, data.frame(class=i, m=m, ms=m+s, m2s=m+2*s, max=max))
}

col<-rainbow(length(levels(group)))

# by median
median <- tapply(number[class=="UCR"], group[class=="UCR"], median)
o <- order(-median, levels(group), decreasing=T)
group <- factor(group, levels=levels(group)[o])
col <- col[o]
#labels <- levels(group)

scale <- function (v) log(1+v)

col2 <- rep(col, each=2)

tiff("boxplotarrangeadtw13hght10nodigitperMB_2cases_3.tiff", width=20*300, height=24*300, res=300, compression="lzw")
#pdf("drawing.pdf", width=8, height=10)

at <- c(0, 10, 100, 600, 1000)

par(mar=c(5,30,1,1)+.1, font=2)

boxplot(scale(number)~class+group, 
		horizontal=T,
		axes=F,
		outpch=16,
		xlab="",
		ylab="",
		#boxcol=col2,
		boxfill=col2,
		#medcol=col2,
		outcol=col2,
		ylim=scale(range(at))
)



labels <- c(t(matrix(c(paste0(levels(group), " (Coding)"),
								paste0(levels(group)," (UCE, n=", table(group[class=="UCR"]), ")",
										"(P",plegend(sapply(levels(group), function (g) p$pvalue[p$type==g])),")")), ncol=2)))



axis(1, at=scale(at), labels=at, font.axis=2, lwd=4, cex.axis=2)
axis(2, at=seq(labels), labels=rep("", length(labels)), las=1, font.axis=2, lwd=4, cex.axis=1.75)

mtext(labels, side=2, line=1, at=seq(labels), las=1, cex=2)

mtext("# UCE mutations per MB per tumor", side=1, line=3, las=1, cex=2.5)

lcol <- c("darkblue", "black")
for (i in 1:nrow(stat)) {
	abline(v=scale(stat[i,2:4]), lty=c("dotted", "dotdash","longdash"), lwd=2, col=lcol[i])
}

stat$class[2]


for (i in 1:nrow(stat)) {
	legend(5, 32-8*i,
			lty=c("dotted", "dotdash","longdash", "blank"),
			lwd=3,
			legend=paste0(c("mean=", "mean+SD=", "mean+2SD=", "max="), 
					formatC(as.numeric(stat[i,2:5]), format="f", digits=0)),
			cex=2,
			bty="n",
			col=lcol[i],
			title.col=lcol[i],
			text.col=lcol[i],
			title=ifelse(i==1, as.character(stat$class[i]), "UCE"),
			xpd=NA)
}

dev.off()


