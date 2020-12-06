


x <- read.table("CLLCRC_SurvAgeGenderNoTreatmPrev.txt", sep="\t", quote="\"", head=T, check.names=F)

number <- x$count

group <- x$DISEASE


m<-mean(number)
s<-sd(number)
max<-max(number)
p <- c(m, m+s, m+2*s)

col<-"black"

group <- factor(group)
col <- rev(col)
labels <- levels(group)

scale <- function (v) log(1+v)

tiff("fig5A.tiff", width=12*300, height=9*300, res=300, compression="lzw")

par(mar=c(3,8,1,9)+.1, font=2)

at <- c(20, 200, 2000)

boxplot(scale(number)~group, 
		horizontal=T,
		axes=F,
		outpch=16,
		boxcol="black",
		boxfill="grey90",
		medcol="black",
		outcol=col,
		outline=T,
		boxwex=.5,
		cex=2.8,
		lwd=4,
		ylim=scale(range(at))
)


axis(1, at=scale(at), labels=at, font.axis=2, lwd=4, mgp=c(3,1.5,0), cex.axis=2.8)
axis(2, at=seq(labels), labels=rep("", length(labels)), las=1, mgp=c(3,1.5,-1.25), font.axis=2, lwd=4, cex.axis=2.8)
mtext(paste0(labels,"\n","(n=", table(group), ")"), side=2, line=3, at=seq(labels), col=col, las=1, cex=2.8, adj=.5)

mtext("# UCR mutations per tumor", side=1, line=5, las=1, cex=2.8)

abline(v=scale(p), lty=c("dotted", "dotdash","longdash"), lwd=3)

legend(5.7, 1.5,
		lty=c("dotted", "dotdash","longdash", "blank"),
		lwd=4,
		legend=paste0(c("mean=", "mean+SD=", "mean+2SD=", "max="), 
				formatC(c(p,max), format="f", digits=0)),
		cex=2.8,
		bty="n",
		xpd=NA)

legend(6,2.5, "p<0.0001",cex=2.8, bty="n")

dev.off()

