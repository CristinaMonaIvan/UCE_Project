

source("source.R")

library(coin)

x <-  read.table("clinicexpression.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F,na="",check.names=F)

time <- "os_months"
event <- "os_status"

x <- x[!is.na(x[,time]),]
x <- x[!is.na(x[,event]),]

x<-x[x[,time]>0, ]

mir <- colnames(x)[20]

d <- NULL

for (m in mir) {
	coxmodel<-coxph(Surv(x[,time], x[,event]) ~ x[,21] + x[,m])
	s <- summary(coxmodel)
	d <- rbind(d, cbind(as.data.frame(
							c("Age",m)),
					as.data.frame(s$conf.int),
					as.data.frame(s$coefficients[,5]),
					t(as.data.frame(s$waldtest))))
					samples=length(which(!is.na(x$Age)))))
	
}

d <- d[order(d$pvalue),]
#d<-subset(d, pvalue<0.05)


write.table(d, file="os_COX_mut_Age.txt", sep="\t", quote=F, row.names=F)


