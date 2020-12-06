

xrel <- function (x) par("usr")[1]+x*diff(par("usr")[1:2])
yrel <- function (y) par("usr")[3]+y*diff(par("usr")[3:4])

comp <- c(A="U", C="G", G="C", T="A", U="A")

jitter.vert2 <- function (v, y0, dist.x, dist.y) {
	s <- NULL
	for (p in v) {
		r <- s[,2][abs(p-s[,1])<dist.x]
		q <- if (length(r)>0) min(setdiff(1:(max(r)+1), r)) else 1
		s <- rbind(s, cbind(p, q))
	}
	y0-dist.y*(s[,2]-1)
}

x1 <- read.table("MIR142completed.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)
x2 <- read.table("bed_filtered_merged_mutect_pindel_568_samples_snpfiltered_ucr_CRC_CLL_UCR5578_noTreatm.txt", sep="\t", head=T, quote="\"", check.names=F, stringsAsFactors=F)

f <- identity

x <- rbind(data.frame(POS=f(x1$POS), REF=x1$REF, ALT=x1$ALT, TYPE=x1$type_pie, TOP=T, stringsAsFactors=F),
		data.frame(POS=f(x2$start), REF=x2$ref_allele, ALT=x2$alt_allele, TYPE=x2$`DISEASE TYPE`, TOP=F, stringsAsFactors=F))

x <- x[!is.na(x$POS),] 
x <- x[order(x$POS),]
x$TYPE <- sub("CLLE", "CLL", x$TYPE)

r <- data.frame(name=c("mir-142", "miR-142-3p", "miR-142-5p", "seed", "seed"),
		start=f(c(56408593, 56408606, 56408644, 56408621, 56408657)),#7meer
		#start=f(c(56408593, 56408606, 56408644, 56408622, 56408658)),#6meer
		end=f(c(56408679, 56408628, 56408664, 56408627, 56408663)),
		stringsAsFactors=F)


col.rect <- c("darkslategray2", "dodgerblue4","deepskyblue3", "red")
names(col.rect) <- unique(r$name)
density.rect <- c(NA, NA, NA, NA)
names(density.rect) <- unique(r$name)

top <- x$TOP
pos <- x$POS
mut <- paste0(x$REF, ">", x$ALT)
alt <- comp[x$ALT]

type <- as.character(x$TYPE)

height <- rep(NA, length(pos))
height[top] <- jitter.vert2(pos[top], .15, 2, -.15)
height[!top] <-jitter.vert2(pos[!top], -.15, 2, +.15)

col.mut <- rainbow(length(unique(mut)))
names(col.mut) <- unique(mut)

pch.type <- 21:22
names(pch.type) <- sort(unique(type))

cex <- 1.5

tiff("MIR142.tiff", width=12*300, height=3*300, res=300, compression="lzw")
par(mar=c(3,2,3,2)+.1,font=2)

plot(NULL, xlim=range(pretty(c(r$start, r$end))), ylim=c(-.4, .95), xlab="", ylab="", cex.main=3, axes=F)
axis(1, at=pretty(c(r$start, r$end)), lwd=3,cex.axis=cex, font.axis=2)

seq <- "GUGUCAUGUGAGUAGGUAUUUCAUCCUUUGUGAUGUGGGAGGUCACGACAAUCAUCACGAAAGAUGAAAUACCCACUGACGUGACAG"
nuc <- strsplit(seq, "")[[1]]
ind <- f(r$start[r$name=="mir-142"])+seq(nuc)-1
col <- rep(NA, length(ind))
for (i in order(r$end-r$start, decreasing=T)) col[r$start[i]<=ind & ind<=r$end[i]] <- col.rect[r$name[i]]
seq.y <- 0
#for (i in seq(nuc)) 
text(ind, seq.y, nuc, cex=.9, col=col, font=2)

segments(x0=pos, y0=ifelse(top, seq.y+.05, seq.y-.05), y1=height, col="darkgrey")

points(pos, height, pch=pch.type[type], bg="white", cex=2.2)

text(pos, height, alt, font=2, cex=.9)

legend(xrel(.25), 1.5, 
		legend=names(col.rect),
		fill=col.rect, density=density.rect, cex=cex, bty="n", xpd=NA, horiz=T)

legend(xrel(.03), 1.5, legend=names(pch.type), pch=pch.type, 
		cex=cex, bty="n", col="black", pt.bg="white", xpd=NA)

legend(xrel(-.05), -.0,
		legend=paste("MDACC Collab"), bty="n", cex=cex, xpd=NA, adj=0)

legend(xrel(-.05), .5, 
		legend=paste("ICGC"), bty="n", cex=cex, xpd=NA, adj=0)

dev.off()

