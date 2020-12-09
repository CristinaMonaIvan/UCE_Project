


library(plotrix)

source("sourcepie3D.R")


x <- read.table("bed_type.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F, na="")

l <- sort(table(x$mut_type))

labels <- names(l)
percentage <- paste0(formatC(100*l/sum(l), format="f", digits=2), "%")

n <- length(l)

cols <- colors()[c(109, 84,508,81,11,494,471, 132)]

labelrad <- rep(1, length(l))
labelrad[labels=="INS"] <- 1.25
labelrad[labels=="DEL"] <- 1.35
labelrad[labels=="C_G"] <- 1.25
labelrad[labels=="A_T"] <- 1.25
labelrad[labels=="A_C"] <- 1.5
labelrad[labels=="A_G"] <- 1.25
labelrad[labels=="C_A"] <- 1.45
labelrad[labels=="C_T"] <- 1.5

labels <- sub("_", ">", labels)

tiff("PieVariantClass2.tif", width=24*300, height=12*300, res=300, compression="lzw")

par(font=2)

ll <- sapply(l, function (n) max(n, sum(l)/100))

pie3D(ll, labels=percentage, col=cols, cex=3, labelcex=2.5, explode=.3, theta=pi/3.1, start=pi*.4, labelrad=labelrad)


legend(1.6, .9, legend=rev(labels), pch=15, col=rev(cols), bty="n", cex=3, xpd=NA)

dev.off()
