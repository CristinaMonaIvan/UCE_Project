



source("mozaicplot.R")


z <- read.table("mergedforFisher.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F, na="", check.names=F)

z$group <- "ChIP-seq-No"
z$group[!is.na(z$chr_hg) & is.na(z$expression)] <- "ChIP-seq-Yes&\nIn vivo-No"
z$group[!is.na(z$chr_hg) & !is.na(z$expression)] <- "ChIP-seq-Yes&\nIn vivo-Yes"

t <- table(z$mut_typeALL, z$group)

col <- c("firebrick1", "lightblue", "green")
names(col) <- colnames(t)

p <- fisher.test(t)$p.value

rownames(t) <- paste0(c("Mut=0", "Mut>0"), "\n(n=", rowSums(t), ")")
colnames(t) <- rep("", ncol(t))

rect.label <- t(apply(t, 1, function (v) formatC(v/sum(v), digits=2, format="f")))

tiff("FisherMutnoMut3grp.tiff", width=2500, height=2700, res=300, compression="lzw")

par(font=2, mar=c(5,15,4,2)+.1, cex.main=2.2)

mosaicplot(t, off=c(1,0), xlab="", ylab="", main="", , cex=2.2, color=col,
		rect.label=rect.label, rect.cex=2.2, rect.font=2)

legend(-.7, .7, legend=names(col), pch=15, col=col, bty="n", xpd=T, cex=2.2,y.intersp=2)


legend(-.65, .95, legend=paste0("p=", formatC(p, digits=2)), xpd=T, bty="n", cex=2.2)

dev.off()
