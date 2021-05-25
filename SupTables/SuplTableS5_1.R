

library(GenomicRanges)

x <- read.table("seq_hg19.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
y <- read.table("coding_intervals_disjoint.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

x <- x[,c("chr_hg19", "start_hg19", "end_hg19")]

colnames(x) <- c("chr", "start", "end")
x$chr <- sub("chr", "", x$chr)

uce <- GRanges(x$chr, IRanges(x$start, x$end))
coding <- GRanges(y$chr, IRanges(y$start, y$end))
uce_coding <- intersect(uce, coding)
uce_noncoding <- setdiff(uce, coding)

r <- rbind(data.frame(type="uce", x),
		data.frame(type="coding", y),
		data.frame(type="uce_coding",
				chr=seqnames(uce_coding),
				start=start(uce_coding),
				end=end(uce_coding)),
		data.frame(type="uce_noncoding",
				chr=seqnames(uce_noncoding),
				start=start(uce_noncoding),
				end=end(uce_noncoding)))

write.table(r, file="regions.txt", sep="\t", quote=F, row.names=F, na="")



