





setwd("C:\\Users\\CAntonescu\\Desktop\\nucleotidechangedsitrib2groups")


cols <- colors()[c(30,109, 81, 508, 84, 11, 471, 494, 132)]


labels<-c("complex","INS","DEL","C>G", "A>T", "A>C", "A>G", "C>A", "C>T")

tiff("legend.tif", width=3*300, height=6*300, res=300, compression="lzw")

plot.new()

legend(.4, .4, xjust=.5, yjust=.5, pch=15, col=rev(cols), legend=rev(labels), pt.cex=5.5, text.font=2, cex=2, bty="n", xpd=NA)

dev.off()




 
		
		