
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

get.region <- function (s, start, name, pattern) {
	g <- gregexpr(pattern, s)
	data.frame(name=name, 
			start=start+g[[1]]-1, 
			end=start+g[[1]]+attr(g[[1]], "match.length")-2)
}

library(Biostrings)
library(RColorBrewer)


x1 <- read.table("MIR99AHG.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
x2 <- read.table("bed_filtered_merged_mutect_pindel_568_samples_201002_snpfiltered_ucr_CRC_CLL_miR99aHG.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

x <- rbind(data.frame(POS=x1$POS, REF=x1$REF, ALT=x1$ALT, TYPE=x1$type_pie, 
				REGION=x1$gene, REGION_START=x1$ucr_start, REGION_END=x1$ucr_end, TOP=T, stringsAsFactors=F),
		data.frame(POS=x2$start, REF=x2$ref_allele, ALT=x2$alt_allele, TYPE=x2$`DISEASE TYPE`,
				REGION=x2$name, REGION_START=x2$start.1, REGION_END=x2$end.1, TOP=F, stringsAsFactors=F))

s <- readBStringSet("UCSCgenomebrowser99AHG.txt")

names(s)
seq.start <- 17442842
seq <- paste(as.character(s), collapse="")

r <- rbind(get.region(seq, seq.start, "intron", "[a-z]+"),
		get.region(seq, seq.start, "exon", "[A-Z]+"),
		get.region(seq, seq.start, "miR-99a-5p", tolower(gsub("U", "T", "AACCCGUAGAUCCGAUCUUGUG"))),
		get.region(seq, seq.start, "miR-99a-3p", tolower(gsub("U", "T", "CAAGCUCGCUUCUAUGGGUCUG"))))

x <- x[order(x$POS),]

x$REF <- as.character(x$REF)
i <- x$REF=="-"
x$REF[i] <- toupper(sapply(x$POS[i], function (j) substr(seq, j-seq.start, j-seq.start)))
x$ALT[i] <- paste0(x$REF[i], x$ALT[i])

u <- x[,c("REGION", "REGION_START", "REGION_END")]
u <- u[!duplicated(u),]
colnames(u) <- c("name", "start", "end")
col.region <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))[1:nrow(u)]
names(col.region) <- u$name

r$name <- as.character(r$name)
col.rect <- c("black", "black", "black", "lightblue")
names(col.rect) <- unique(r$name)[c(4,3,2,1)]
pch.rect <- c(8,3,4,15)
names(pch.rect) <- names(col.rect)

top <- x$TOP
pos <- x$POS
mut <- paste0(x$REF, ">", x$ALT)

type <- as.character(x$TYPE)
region <- as.character(x$REGION)

height <- rep(NA, length(pos))
height[top] <- jitter.vert2(pos[top], 1.2, 10000, .15)
height[!top] <-jitter.vert2(pos[!top], -1.5, 10000, -.15)

col.mut <- rainbow(length(unique(mut)))
names(col.mut) <- unique(mut)

pch.type <- c(12:20, 23, 25)[1:length(unique(type))]
names(pch.type) <- sort(unique(type))

tiff("MIR99AHG.tiff", width=12*300, height=6*300, res=300, compression="lzw")
par(mar=c(4,3,6,8)+.1,font=2)

plot(NULL, xlim=range(pretty(c(r$start, r$end))), ylim=c(-1.5, 1.5), xlab="", ylab="", cex.main=3, axes=F)
axis(1, at=pretty(c(r$start, r$end)), lwd=3, cex.axis=1.5)

rect.lower <- -.1
rect.upper <- .1

segments(x0=pos, y0=ifelse(top, rect.upper, rect.lower), y1=height, col=col.region[region], lwd=2)

rect(r$start-1, rect.lower, r$end, rect.upper, col=col.rect[r$name], border=col.rect[r$name])

points(pos, height,
		col=col.mut[mut], bg=col.mut[mut], pch=pch.type[type], cex=2.2)

sym <- r[r$name!="intron"&!r$name%in%u$name,]
sym$mid <- (sym$start+sym$end)/2
sym <- sym[order(sym$mid),]

sym.height <- jitter.vert2(sym$mid, rect.upper-.01, 2000, -.1)
points(sym$mid, sym.height, pch=pch.rect[sym$name], cex=1.5, lwd=2)

legend(xrel(-.02), 2.55, 
		legend=names(col.mut),
		pch=16, cex=1.2, bty="n", col=col.mut, xpd=NA)

legend(xrel(.18), 2.55, legend=names(pch.type), pch=pch.type, 
		cex=1.2, bty="n", col="black", pt.bg="black", xpd=NA)

legend(xrel(.37), 2.55, 
		legend=names(col.rect),
		pch=pch.rect, cex=1.2, bty="n", col=col.rect, xpd=NA)

legend(xrel(.6), 2.55, 
		legend=names(col.region),
		pch=15, cex=1.2, bty="n", col=col.region, xpd=NA, ncol=3)

legend(xrel(1.08), -.4,
		legend=paste("MDACC Collab"), bty="n", cex=1.5, xpd=NA, adj=1)

legend(xrel(1.08), .75, 
		legend=paste("ICGC"), bty="n", cex=1.5, xpd=NA, adj=1)

dev.off()



