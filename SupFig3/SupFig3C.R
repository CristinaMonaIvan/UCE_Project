



library(barplot3d)
library(rgl)

x <- read.table("mut_perGene_r.txt", sep="\t", head=T, quote="\"", fill=T,row.names=1)

types <- grep("percentage", colnames(x))[1:22]

M <- as.matrix(x[order(x$total_percentage, x$max_percentage, decreasing=T)[1:30],types])

write.table(M, file="data.txt", sep="\t", quote=F, row.names=T, na="")

col <- rep(rainbow(22), 30)
col[t(M)==0] <- "white"

colnames(M) <- sub("_percentage", "", colnames(M))

rownames(M)[3] <- sub(".*\\|", "", rownames(M)[3])

par(mar=c(5, 4, 4, 2) + 0.1, xpd=NA)



barplot3d(rows=30,cols=22,z=t(M),scalexy=5, gap=0.3, alpha=0.2,theta=45,phi=50,
topcolors=col, xlabels = colnames(M),ylabels=rownames(M),
xsub="",ysub="",zsub="")

rgl.snapshot("plot.png", fmt = "png", top = TRUE)





