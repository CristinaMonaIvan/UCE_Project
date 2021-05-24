



library(RColorBrewer)

library(alluvial)

source("C:\\Users\\CAntonescu\\eclipse-workspace\\ALLUVIALDiagrams\\sourcepertype.R")

chrorder <- function (chr) as.numeric(sub("chr", "", sub("X", 23, sub("Y", 24, chr))))

x <- read.table("tableuceinfo.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
w <- read.table("indelsnv_seq_hg19_allcoding_perpat.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

uce <- read.table("uce.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
uce$UCE <- paste0("UCE_", uce$UCE)
uce$label <- paste0(uce$UCE, " (", uce$gene, ")")

w$type[w$type=="PAAD"] <- "PACA"
w$type[w$type=="SARC"] <- "BSTC"
w$type[w$type=="Glioma"] <- "CNSC"
w$type[w$type=="CLLE"] <- "CLL"

typecol <- rainbow(length(unique(x$type)))
names(typecol) <- sort(unique(x$type))

#t <- tapply(x$total, x$type, sum)
t <- tapply(w$count_UCR, w$type, median)

x <- x[order(-t[x$type], x$type, chrorder(x$chr_hg19), x$start_hg19, decreasing=T),]

for (type in unique(x$type)) {
	
	tiff(paste0(type, ".tif"), width=4000, height=8000, res=300, compression="lzw")
	
	a <- alluvial(x[,c("type", "chr_hg19")], freq=x$total,
			ordering=list(order(-t[x$type], x$type, chrorder(x$chr_hg19), x$start_hg19, decreasing=T),
					order(chrorder(x$chr_hg19), x$start_hg19, decreasing=T)),
			gap.width=0, axis_labels=rep("", 2), cex=c(2, 2),
			cw = c(.12, .08),
			col=typecol[x$type],
			rect.col=list(typecol[x$type], NULL),
			alpha=1, border=NA, labels=T,
			hide=x$type!=type,
			show.blocks = list(type, NULL),
			mar=c(2, 1, 1, 12))
	
	v <- data.frame(UCE=x$UCE, y=apply(a$endpoints[[2]], 1, mean))
	v <- v[v$UCE%in%uce$UCE,]
	v <- aggregate(y~UCE, v, mean)
	v <- merge(v, uce)
	
	text(2.1, v$y, v$label, xpd=NA, adj=c(0, .5))
	segments(2.08, v$y, 2.09, v$y, xpd=NA)
	
	dev.off()
}

