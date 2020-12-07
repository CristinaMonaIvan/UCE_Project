



source("mozaicplot.R")


z <- read.table("clinicsurvCompletemut.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F, na="")

tiff("ResponseLowHyper.tiff", width=2500, height=2700, res=300, compression="lzw")

par(font=2, mar=c(5,9.5,4,2)+.1, cex.main=2.2)

z$type_mut_2<-factor(z$type_mut_2)


levels(z$type_mut_2) <- c("lowly&\nnon","highly&\nhyper")

z$type_mut_2 <- factor(z$type_mut_2, levels=c("highly&\nhyper", "lowly&\nnon"))


z <- z[z$disease_status_last_followup!="",]

z$disease_status_last_followup<-factor(z$disease_status_last_followup)
#x$disease_status_last_followup[x$disease_status_last_followup=="partial remission"] <- "stable"

levels(z$disease_status_last_followup)<-c("CR","NED","PR","PD","Relap","SD")

z$disease_status_last_followup <- factor(z$disease_status_last_followup, levels=c("NED", "CR", "PR", "SD", "Relap", "PD"))

fisher.test(table(z$type_mut_2, z$disease_status_last_followup), simulate.p.value=T)

fn<- table(z$type_mut_2, z$disease_status_last_followup)
fp<-fisher.test(fn, simulate.p.value=T)$p.value
y<- data.frame(pvalue=fp, t(c(fn)))


t <- table(z$type_mut_2, z$disease_status_last_followup)

f <- colnames(t)
colnames(t) <- rep("", ncol(t))

col <- rainbow(ncol(t))

rect.label <- t(apply(t, 1, function (v) formatC(v/sum(v), digits=2, format="f")))

mosaicplot(t, off=c(1,0), xlab="", ylab="", main="", , cex=2.2, color=col,
		rect.label=rect.label, rect.cex=2.2, rect.font=2)

legend(-.3, .6, legend=f, pch=15, col=col, bty="n", xpd=T, cex=2.2)


legend(-.45, .95, legendplegend(fisher.test(t, simulate.p.value=T)$p.value), xpd=T, bty="n", cex=2.2)

dev.off()

