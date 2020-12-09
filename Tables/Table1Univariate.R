



source("source.R")

library(coin)



x <-  read.table("clinicexpression.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F,na="",check.names=F)


time <- "os_months"
event <- "os_status"

x <- x[!is.na(x[,time]),]
x <- x[!is.na(x[,event]),]

x<-x[x[,time]>0,]


mir <- colnames(x)[c(8,10,11,20)]


d <- NULL

for (m in mir) {
	y<-x[!is.na(x[,m]),]
	coxmodel<-coxph(Surv(y[,time], y[,event])~ y[,m])
	s <- summary(coxmodel)
	z <- cox.zph(coxmodel)$table
	d <- rbind(d, cbind(data.frame(mir=m), as.data.frame(s$conf.int), t(as.data.frame(s$waldtest)), pt=length(which(!is.na(x[,m]))), t(z[1,-2])))
	}

d$fdr <- p.adjust(d$pvalue, method="fdr")
#d<-subset(d, pvalue<0.05)
#d<-d[d$"exp(coef)">1,]



write.table(d, file="os_COX_univariate.txt", sep="\t", quote=F, row.names=F)

