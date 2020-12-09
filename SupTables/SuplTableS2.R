

tables <- function (w) {
	d <- NULL
	for (i in colnames(w)) {
		if (class(w[,i])=="factor" || length(unique(w[,i]))<=3) {
			t <- table(w[,i], useNA="ifany")
			d <- rbind(d, data.frame(var=c(i, names(t), ""), 
							number=c("", paste0(t, " (", formatC(100*t/sum(t), format="f", digits=1), ")"), "")))
		}
		else
		{
			d <- rbind(d, data.frame(var=c(i, "range", "median", ""),
							number=c("", 
									paste0("(", min(w[,i], na.rm=T), ",", max(w[,i], na.rm=T), ")"),
									median(w[,i], na.rm=T),
									"")))
		}
	}
	d
}


x <- read.table("clinicsurvmut.txt", sep="\t", head=T, na="")

s <- tables(x[,3:5])

write.table(s, file="GenderAge.txt", sep="\t", row.names=F, col.names=F, quote=F, na="NA")

for (t in unique(x$type)) {
	s <- tables(x[x$type==t,3:5])
	write.table(s, file=paste0("GenderAge_", t, ".txt"), sep="\t", row.names=F, col.names=F, quote=F, na="NA")
}
