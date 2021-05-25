

x <- read.table("mutevents.txt", sep="\t", head=T, quote="\"", check.names=F, na="",stringsAsFactors=F)
y <- read.table("mut_regions.txt", sep="\t", head=T, quote="\"", check.names=F, na="",stringsAsFactors=F)

#x$REF_length <- nchar(x$REF)
#x$ALT_length <- nchar(x$ALT)

z <- merge(x, y, sort=F)

write.table(z, file="overlap.txt", sep="\t", quote=F, row.names=F)

for (r in unique(z$type)) {
	write.table(z[z$type==r,colnames(z)!="type"], file=paste0("overlap_", r, ".txt"), sep="\t", quote=F, row.names=F)
}