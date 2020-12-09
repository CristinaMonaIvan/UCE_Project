

xrel <- function (x) par("usr")[1]+x*diff(par("usr")[1:2])
yrel <- function (y) par("usr")[3]+y*diff(par("usr")[3:4])

jitter.vert2 <- function (v, y0, dist.x, dist.y) {
	s <- NULL
	for (p in v) {
		r <- s[,2][abs(p-s[,1])<dist.x]
		q <- if (length(r)>0) min(setdiff(1:(max(r)+1), r)) else 1
		s <- rbind(s, cbind(p, q))
	}
	y0-dist.y*(s[,2]-1)
}

x1 <- read.table("ICGC_3075.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
x2 <- read.table("bed_filtered_merged_mutect_pindel_568_samples_snpfiltered_ucr_CRC_CLL_UCR3075_noTreatm.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

r <- data.frame(name=c("CNOT2 (NM_014515) 3'UTR", "UCE_3075"),
		start=c(70747696, 70748292),
		end=c(70748773, 70748630))

col.rect <- c("darkslategray3", "deepskyblue4")

x0 <- x1[x1$gene=="UCR_3075",]
x0 <- x0[order(x0$POS),]

pos <- x0$POS
mut <- paste0(x0$REF, ">", x0$ALT)

type <- factor(x0$type_pie)
type <- as.character(type)

height <- jitter.vert2(pos, 1.0, 20, .12)

col.mut <- rainbow(8)[1:7]
names(col.mut) <- unique(mut)

pch.type <- c(15:18,20,23,25)
names(pch.type) <- sort(unique(type))

tiff("UCR_3075_all_noTreat.tiff", width=12*300, height=5*300, res=300, compression="lzw")
par(mar=c(5,5,4,5)+.1,font=2)

plot(NULL, xlim=range(pretty(c(r$start, r$end))), ylim=c(-.7, 1.5), xlab="", ylab="", cex.main=3, axes=F)
axis(1, at=pretty(c(r$start, r$end)),lwd=3,cex.axis=1.7, font.axis=2)

rect.lower <- -.1
rect.upper <- .1

rect(r$start-1, rect.lower, r$end, rect.upper, col=col.rect)

segments(x0=pos, y0=rect.upper, y1=height, col="darkgrey", lwd=2)

points(pos, height,
		col=col.mut[mut], bg=col.mut[mut], pch=pch.type[type], cex=2)

legend(xrel(-.05), 1.95, 
		legend=names(col.mut),
		pch=16, cex=1.5, bty="n", col=col.mut, xpd=NA)


legend(xrel(.3), 1.95, legend=names(pch.type), pch=pch.type, 
		cex=1.5, bty="n", col="black", pt.bg="black", xpd=NA)

legend(70748300, 1.95, 
		legend=r$name,
		pch=15, cex=1.5, bty="n", col=col.rect, xpd=NA)

# MDACC

x0 <- x2

pos <- x0$start
mut <- paste0(x0$ref_allele, ">", x0$alt_allele)

type <- factor(x0$`DISEASE TYPE`)
type <- as.character(type)

height <- jitter.vert2(pos, -.5, 20, .12)

col.mut <- rainbow(8)[8]
names(col.mut) <- unique(mut)

pch.type <- c(17)
names(pch.type) <- sort(unique(type))

segments(x0=pos, y0=rect.lower, y1=height, col="darkgrey", lwd=2)

points(pos, height,
		col=col.mut[mut], bg=col.mut[mut], pch=c(25), cex=2)


legend(xrel(-.05), -.1, 
		legend=names(col.mut),
		pch=16, cex=1.5, bty="n", col=col.mut, xpd=NA)

legend(xrel(.3), -.1, 
		legend=names(pch.type),
		pch=pch.type, 
		cex=1.5, bty="n", col="black", pt.bg="black", xpd=NA)

legend(xrel(1), -.27, 
		legend=paste("MDACC Collab"), bty="n", cex=1.5, xpd=NA, adj=1)


legend(xrel(1), .8, 
		legend=paste("ICGC"), bty="n", cex=1.5, xpd=NA, adj=1)

dev.off()

