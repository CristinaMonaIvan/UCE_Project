


library(plotrix)

source("sourcepie3D.R")


x <- read.table("allvariantslow.txt", sep="\t", head=T, quote="\"", stringsAsFactors=F, na="")

l <- sort(table(x$change))

labels <- names(l)
percentage <- paste0(formatC(100*l/sum(l), format="f", digits=2), "%")

n <- length(l)

cols <- colors()[c(30,109, 81, 11,508,84,494,471,132)]


labelrad <- rep(1, length(l))
labelrad[labels=="complex"] <- 1.2
labelrad[labels=="INS"] <- 1.25
labelrad[labels=="DEL"] <- 1.15
labelrad[labels=="CtoG"] <- 1.5
labelrad[labels=="AtoT"] <- 1.5
labelrad[labels=="AtoC"] <- 1.5
labelrad[labels=="AtoG"] <- 1.35
labelrad[labels=="CtoA"] <- 1.45
labelrad[labels=="CtoT"] <- 1.5

labels <- sub("to", ">", labels)

tiff("PieVariantClassLow.tif", width=24*300, height=12*300, res=300, compression="lzw")

par(font=2)

ll <- sapply(l, function (n) max(n, sum(l)/100))

pie3D(ll, labels=percentage, col=cols, cex=3, labelcex=2.5, explode=.3, theta=pi/3.1, start=pi*.4, labelrad=labelrad)


dev.off()

