

#setwd("C:\\Users\\CAntonescu\\Desktop\\workforUCRpaper\\UCRegions\\UCR_events\\types\\perpat")

setwd("C:\\Users\\CAntonescu\\Desktop\\try")


x <- read.table("types.txt", sep="\t", quote="\"", head=T, check.names=F)


number <- c(x$count_lincRNAandRNA,
		x$count_Intron,
		x$count_IGR,
		x$count_CDS,
		x$count_3pand5p)

group <- factor(rep(x$type, 5))
class <- rep(1:5, each=nrow(x))
classnames <- c("lincRNA+RNA", "Intron", "IGR", "CDS", "3p+5p")

stat <- NULL
for (i in unique(class)) {
	m<-mean(number[class==i])
	s<-sd(number[class==i])
	max<-max(number[class==i])
	stat <- rbind(stat, data.frame(class=i, m=m, ms=m+s, m2s=m+2*s, max=max))
}

col<-rainbow(length(levels(group)))

# by median
median <- tapply(number, group, median)
o <- order(-median, levels(group), decreasing=T)
group <- factor(group, levels=levels(group)[o])
col <- col[o]
#labels <- levels(group)

scale <- function (v) log(1+v)

col2 <- rep(col, each=5)

tiff("VariantClassification_perCancerType_2.tiff", width=21*300, height=45*300, res=300, compression="lzw")
#pdf("drawing.pdf", width=8, height=10)

at <- c(0, 10, 100, 1000,2000)

par(mar=c(5,24,1,2)+.1, font=2)

boxplot(scale(number)~class+group, 
		horizontal=T,
		axes=F,
		outpch=16,
		xlab="",
		ylab="",
		boxfill=col2,
		outcol=col2,
		cex=2,
		lwd=2,
		ylim=scale(range(at))
)


labels <- c(t(matrix(c(paste0(levels(group), " (lincRNA+RNA)"),
								paste0(levels(group), " (Intron)"),
								paste0(levels(group), " (IGR)"),
								paste0(levels(group), " (CDS)"),
								paste0(levels(group)," (3p+5p, n=", table(group[class==1]), ")")), ncol=5)))

axis(1, at=scale(at), labels=at, font.axis=2, lwd=6, cex.axis=4, mgp=c(3,3,0))
axis(2, at=seq(labels), labels=rep("", length(labels)), las=1, font.axis=2, lwd=6, cex.axis=2)

mtext(labels, side=2, line=1, at=seq(labels), font=2, las=1, cex=2)

#mtext("# UCR mutations per tumor", side=1, font=2, line=3, las=1, cex=2)

lcol <- c("darkblue", "darkgreen", "black", "magenta", "darkred")
for (i in 1:nrow(stat)) {
	abline(v=scale(stat[i,2:4]), lty=c("dotted", "dotdash","longdash"), lwd=2, col=lcol[i])
}

for (i in 1:nrow(stat)) {
	legend(5, 9*i,
			lty=c("dotted", "dotdash","longdash", "blank"),
			lwd=3,
			legend=paste0(c("mean=", "mean+SD=", "mean+2SD=", "max="), 
					formatC(as.numeric(stat[i,2:5]), format="f", digits=0)),
			cex=3,
			bty="n",
			col=lcol[i],
			title.col=lcol[i],
			text.col=lcol[i],
			title=classnames[stat$class[i]],
			xpd=NA)
}

dev.off()


