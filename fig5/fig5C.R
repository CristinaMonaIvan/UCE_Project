


library(gplots)

x <- read.table("mirnaevents.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
w <- read.table("table_type.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

rownames(w) <- w$type

z <- x[,c("gene_symbol", "chr", "gene_start")]
z <- z[!duplicated(z),]
z <- z[order(as.numeric(sub("X", "23", sub("chr", "", z$chr))), z$gene_start),]

z <- z[!duplicated(z$gene_symbol),]

y <- as.matrix(table(x$type_pie, x$gene_symbol))
y <- y[,z$gene_symbol]

chr <- tapply(1:nrow(z), z$chr, range)
chr <- chr[unique(z$chr)]
chr.mid <- sapply(chr, mean)
chr.sep <- sapply(chr, "[", 1)[-1]-.5

tiff("fig5C.tiff", width=3000, height=3900, res=300, compression="lzw")

layout(matrix(1:2), heights=c(10,1))

par(mar=c(10,5,4,9)+.1, cex.axis=1.5)


n <- 5

col <- c("white", colorpanel(n, "green", "purple"))

image(1:nrow(y), 1:ncol(y), pmin(y,n), col=col, zlim=c(0,n), axes=F, xlab="", ylab="")

axis(1, at=1:nrow(y), labels=paste0(rownames(y), "(", w[rownames(y),"Patients"], ")"), font.axis=2,las=2)

odd <- seq(from=1, to=length(chr), by=2)
even <- seq(from=2, to=length(chr), by=2)

axis(2, at=chr.mid, labels=names(chr.mid), las=1, font.axis=2, lwd.ticks=0, mgp=c(3,.5,0))


axis(4, at=1:ncol(y), labels=colnames(y), las=1, lwd.ticks=0, font.axis=4, mgp=c(3,.5,0))


box()

abline(v=2:(ncol(y)-1)-.5, lwd=1.5,lty=2)
abline(h=chr.sep, lwd=1.5,lty=2)

par(mar=c(3,4,1,4))
image(x=0:n, y=1, z=matrix(0:n), col=col, axes=F, xlab="", ylab="")
axis(1, at=0:n, font.axis=2)
box()

dev.off()