


source("source.R")

library(coin)

x <- read.table("clinicsurvCompletemut.txt", sep="\t", head=T, quote="", check.names=F, stringsAsFactors=F,na="")


time <- "os_months"
event <- "os_status"

x <- x[!is.na(x[,time]),]


x <- x[!is.na(x[,event]),]


x<-x[x[,time]>0, ]


x$levels <-factor(x$type_mut)


d <- sapply(c(0,53,106,159,212,265), function (t) table(x$levels[x[,time]>=t]))

colnames(d) <-c(0,53,106,159,212,265)



survtest(x[,time], x[,event], x$levels)

z <- logrank_test(Surv(x[,time], x[,event]) ~ x$levels)


survfit(Surv(x[,time], x[, event]) ~ x$levels)

formatC(summary(survfit(Surv(x[,time], x[, event]) ~ x$levels))$table[,"median"], format="f", digits=1)
tiff("OSlowvshigh.tif", width=2700, height=2700, res=300, compression="lzw")
par(font=2)
par(oma=c(0,2,0,0), mar=c(14,6,4,6)+.1)
plot(survfit(Surv(x[,time], x[, event]) ~ x$levels), xlim=c(0,265), xlab="", ylab="", axes=F, main="", lwd=4, lty=1, 
		cex.lab=2.5, font.lab=2, mark=3, cex=2.5, cex.main=2.5, col=c("red", "blue"))
axis(2, at=c(0, .2,.4, .6, .8, 1), labels=c(0, 20, 40, 60, 80, 100), las=1, lwd=4, cex.axis=2.5, font=2)
mtext("OS (months)", side=1, line=4.5, cex=2.5)
axis(1, at=c(0,53,106,159,212,265), lwd=4, font=2, cex.axis=2.5, mgp=c(3,1.5,0))
mtext("Percentage survival (%)", side=2, line=5, cex=2.5)
text(80, .7, "(77.5)", cex=2.5, col="blue")
text(80, .2, "(30.6)", cex=2.5, col="red")
legend("topright", plegend(survtest(x[,time], x[,event], x$levels)), cex=2.5, bty="n")
mtext(c(d[1,-ncol(d)], ""), at=c(0,53,106,159,212,270), side=1, line=11.5, cex=2.5, col="red")
mtext("highly&\nhyper", at=270, side=1, line=13, cex=2.5, col="red")
mtext(c(d[2,-ncol(d)], "lowly&non"), at=c(0,53,106,159,212,270), side=1, line=8.5, cex=2.5, col="blue")
mtext("# at risk", at=265, side=1, line=4.5 , cex=2.5)
dev.off()

