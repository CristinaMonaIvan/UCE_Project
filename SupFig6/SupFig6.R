

source("source.R")


library(coin)

x <- read.table("mutcodinggrupeclinic.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F)



time <- "os_months"

event <- "os_status"


x$levels <-factor(x$type_codmut_2)


d <- sapply(c(0, 53, 106, 159, 212, 265), function (t) table(x$levels[x[,time]>=t]))

colnames(d) <-c(0, 53, 106, 159, 212, 265)



survtest(x[,time], x[,event], x$levels)

z <- logrank_test(Surv(x[,time], x[,event]) ~ x$levels)


survfit(Surv(x[,time], x[, event]) ~ x$levels)

formatC(summary(survfit(Surv(x[,time], x[, event]) ~ x$levels))$table[,"median"], format="f", digits=1)
tiff("SurvCoding2Groups.tif", width=2700, height=2700, res=300, compression="lzw")
par(font=2)
par(oma=c(0,2,0,0), mar=c(14,6,4,6)+.1)
plot(survfit(Surv(x[,time], x[, event]) ~ x$levels), xlim=c(0,265), xlab="", ylab="", axes=F, main="", lwd=4, lty=1, 
		cex.lab=2.5, font.lab=2, mark=3, cex=2.5, cex.main=2.5, col=c("blue" ,"red"))
axis(2, at=c(0, .2,.4, .6, .8, 1), labels=c(0, 20, 40, 60, 80, 100), las=1, lwd=4, cex.axis=2.5, font=2)
mtext("Percentage survival (%)", side=2, line=5, cex=2.5)
mtext("OS (months)", side=1, line=4.5, cex=2.5)
axis(1, at=c(0, 53, 106, 159, 212, 265), lwd=4, font=2, cex.axis=2.5, mgp=c(3,1.5,0))
legend("topright", "p=plegend(survtest(x[,time], x[,event], x$levels))", cex=2.5, bty="n")
mtext(c(d[1,-ncol(d)], "low&null"), at=c(0, 53, 106, 159, 212, 265), side=1, line=8.5, cex=2.3, col="blue")
mtext(c(d[2,-ncol(d)], "high&hyper"), at=c(0, 53, 106, 159, 212, 265), side=1, line=11.5, cex=2.3, col="red")
mtext("# at risk", at=265, side=1, line=4.5 , cex=2.5)
dev.off()

