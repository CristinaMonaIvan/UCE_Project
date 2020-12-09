

x <- read.table("bed_filtered_merged_mutect_pindel_568_samples_201002_snpfiltered_ucr_CRC_CLL.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

x$group <- ceiling(10*(x$start-x$start_hg19+1)/(x$end_hg19-x$start_hg19+1))

m <- table(x$group)
p <- m/sum(m)*100

tiff("number_MDACC.tiff", height=7*300, width=7*300, res=300, compression="lzw")
b <- barplot(m, axes=F, names=F, ylim=range(pretty(m)))
axis(1, at=seq(from=1.5*b[1]-.5*b[2], by=b[2]-b[1], length.out=length(b)+1),
		labels=paste0((-1):10, "/10"), cex.axis=.8)
axis(2, las=1, at=pretty(range(m)))
dev.off()

tiff("percent_MDACC.tiff", height=7*300, width=7*300, res=300, compression="lzw")
b <- barplot(p, axes=F, names=F, ylim=range(pretty(p)))
axis(1, at=seq(from=1.5*b[1]-.5*b[2], by=b[2]-b[1], length.out=length(b)+1),
		labels=paste0((-1):10, "/10"), cex.axis=.8)
axis(2, las=1, at=pretty(range(p)), labels=paste0(pretty(range(p)), "%"))
dev.off()

d <- data.frame(interval=paste0((-1):9, "/10 - ", 0:10, "/10"),
		count=c(m), percentage=c(p))

write.table(d, file="table_MDACC.txt", sep="\t", quote=F, row.names=F)
