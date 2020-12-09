


xrel <- function (x) par("usr")[1]+x*diff(par("usr")[1:2])
yrel <- function (y) par("usr")[3]+y*diff(par("usr")[3:4])

complement <- function (x) {
	comp <- c(A="T", C="G", G="C", T="A")
	sapply(strsplit(x, NULL), function (y) paste(comp[y], collapse=""))
}

jitter.vert2 <- function (v, y0, dist.x, dist.y) {
	s <- NULL
	for (p in v) {
		r <- s[,2][abs(p-s[,1])<dist.x]
		q <- if (length(r)>0) min(setdiff(1:(max(r)+1), r)) else 1
		s <- rbind(s, cbind(p, q))
	}
	y0-dist.y*(s[,2]-1)
}

region <- function (s, start, name, pattern) {
	g <- gregexpr(pattern, s)
	data.frame(name=name, 
			start=start+g[[1]]-1, 
			end=start+g[[1]]+attr(g[[1]], "match.length")-2)
}

library("Biostrings")

x1 <- read.table("UCR_2905events_ICGC.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
x2 <- read.table("UCR_2905events_MDACC.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

s <- readBStringSet("KRAS.txt")

names(s)
seq.start <- 25358180
seq <- reverse(paste(as.character(s), collapse=""))

r <- rbind(region(seq, seq.start, "KRAS intron", "[a-z]+"),
		region(seq, seq.start, "KRAS exon", "[A-Z]+"))

r <- rbind(r[c(5,10),],
		data.frame(name="UCR 2905", start=25398239, end=25398354))

r$end[1] <- 25398400

r$name <- as.character(r$name)
col.rect <- c("UCR 2905"="black", "KRAS exon"="lightblue", "KRAS intron"="lightgreen")
density.rect <- c("UCR 2905"=20, "KRAS exon"=NA, "KRAS intron"=NA)

x0 <- x1

x0 <- x0[,c("POS", "REF", "ALT", "type")]
x0 <- aggregate(1:nrow(x0), x0, length)
colnames(x0)[5] <- "count"
x0 <- x0[order(x0$POS),]


pos <- x0$POS
mut <- paste0(complement(x0$REF), ">", complement(x0$ALT))

height <- jitter.vert2(pos, 2.5, 3, .1)

type <- x0$type

col.mut <- rainbow(length(unique(mut)))
names(col.mut) <- sort(unique(mut))

pch.type <- 9:20
names(pch.type) <- sort(unique(type))

tiff("UCR2905bothcohorts2.tiff", width=12*300, height=7.5*300, res=300, compression="lzw")
par(mar=c(5,5,4,1)+.1,font=2)

plot(NULL, xlim=range(pretty(c(r$start, r$end))), ylim=c(-.7, 2.5), xlab="", ylab="", cex.main=3, axes=F)
axis(1, at=pretty(c(r$start, r$end)),lwd=3,cex.axis=1.7)

rect.lower <- -.1
rect.upper <- .1

rect(r$start-1, rect.lower, r$end, rect.upper, 
		col=col.rect[r$name], border=col.rect[r$name], density=density.rect[r$name])

segments(x0=pos, y0=rect.upper, y1=height, col="darkgrey",lwd=2)

points(pos, height,
		col=col.mut[mut], bg=col.mut[mut], pch=pch.type[type], cex=2.2)


legend(xrel(-.05), 2.65, 
		legend=paste0(names(col.mut), " (", sapply(names(col.mut), function (m) sum(x0$count[mut==m])), ")"),
		pch=16, cex=1.2, bty="n", col=col.mut, xpd=NA)

legend(xrel(.1), 2.65, 
		legend=paste0(names(pch.type), " (", sapply(names(pch.type), function (m) sum(x0$count[type==m])), ")"),
		pch=pch.type, 
		cex=1.2, bty="n", col="black", pt.bg="black", xpd=NA)

legend(xrel(.72), 2.65, 
		#legend=names(col.rect),
		legend=c("UCE_2905", "KRAS 2nd Exon", "KRAS Intron"),
		density=density.rect,
		cex=1.2, bty="n", fill=col.rect, xpd=NA)

# MDACC

x0 <- x2

x0 <- x0[,c("POS", "REF", "ALT", "TISSUE SITE")]
x0 <- aggregate(1:nrow(x0), x0, length)
colnames(x0)[5] <- "count"
x0 <- x0[order(x0$POS),]


pos <- x0$POS
mut <- paste0(complement(x0$REF), ">", complement(x0$ALT))

height <- jitter.vert2(pos, -.7, 3, -.12)

type <- sub("Colon", "CRC", x0$`TISSUE SITE`)

col.mut <- rainbow(length(unique(mut)))
names(col.mut) <- sort(unique(mut))

pch.type <- c(8,12)
names(pch.type) <- sort(unique(type))

segments(x0=pos, y0=rect.lower, y1=height, col="darkgrey",lwd=2)

points(pos, height,
		col=col.mut[mut], bg=col.mut[mut], pch=pch.type[type], cex=2.2)


legend(xrel(-.05), -.1, 
		legend=paste0(names(col.mut), " (", sapply(names(col.mut), function (m) sum(x0$count[mut==m])), ")"),
		pch=16, cex=1.2, bty="n", col=col.mut, xpd=NA)

legend(xrel(.1), -.1, 
		legend=paste0(names(pch.type), " (", sapply(names(pch.type), function (m) sum(x0$count[type==m])), ")"),
		pch=pch.type, 
		cex=1.2, bty="n", col="black", pt.bg="black", xpd=NA)

legend(xrel(.7), -.2, 
		legend=paste("MDACC Collab"), bty="n", cex=2, xpd=NA)


legend(xrel(.7), .9, 
		legend=paste("ICGC"), bty="n", cex=2, xpd=NA)

dev.off()


