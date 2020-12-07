
pointplot <- function(l, jitter=.15, normal=T, ylab.adj=0, top.adj=0, ylab="", bold=F, bound=Inf, cex.axis=1, group.names=names(l), main="", arrow.col="black", arrow.lwd=2.5, pch=18:16, cex=1.4,...) {
	l <- lapply(l, function (y) y[!is.na(y)])
	u <- lapply(l, function (y) y[y>bound])
	v <- lapply(l, function (y) y[y<=bound])
	r <- pretty(c(do.call(c, v), if (normal) sapply(v, function (y) mean(y)+c(-1, 1)*sd(y)) else NULL))
	par(mar=c(5,5+ylab.adj,4+top.adj+max(sapply(u, length)),2)+.1)
	f <- if (bold) 2 else 1
	par(font=f)
	stripchart(v, vert=T, method="jitter", jitter, pch=pch, ylim=c(min(r), max(r)+.1*diff(range(r))), frame.plot=F, yaxt="n", lwd=f, font.axis=f, cex.axis=cex.axis, cex=cex, group.names=group.names, ...)
	title(main, line=par("mar")[3]-2, ...)
	#axis(2, las=1, at=r, lwd=f, font.axis=f, cex.axis=cex.axis)
	mtext(ylab, side=2, line=4+ylab.adj, cex=cex.axis)
	#mtext(sapply(l, function (y) paste("n=", length(y), sep="")), side=1, line=2.5, at=seq(l), cex=cex.axis)
	sp <- diff(par("usr")[3:4])/20
	for (i in seq(l)) {
		lines(cbind(i+c(-1, 1)*jitter, (if (normal) mean else median)(l[[i]])), lty=2, col=arrow.col)
		if (normal) arrows(x0=i, y0=mean(l[[i]])-sd(l[[i]]), y1=mean(l[[i]])+sd(l[[i]]), angle=90, code=3, length=.1, col=arrow.col, lwd=arrow.lwd)
		points(x=rep(i-.075*length(l), length(u[[i]])), y=par("usr")[4]+(-1.5+seq(as.list(u[[i]])))*sp, xpd=T, pch=18, cex=cex,...)
		if (length(u[[i]])>0) text(x=rep(i, length(u[[i]])), y=par("usr")[4]+(-1.5+seq(as.list(u[[i]])))*sp, labels=formatC(sort(u[[i]]), format="f", digits=3), xpd=T, cex=cex.axis)
	}
}

#jitter=.15
survplot <- function (time, event, group, names=levels(group), footnote=NULL, legend.pos="bottomleft", add.median=F, col=par("col"), legend.space=1, bold=F, add.legend=T, xlab="survival months", ylab="probability of survival", lty=1:4, cex.axis=1, cex.legend=1, cex.mark=1, add.surv50=F, ...) {
	r <- pretty(time)
	if (!is.null(footnote)) par(oma=c(2,0,0,0))
	f <- if (bold) 2 else 1
	par(font=f)
	plot(survfit(Surv(time, event) ~ group), xlim=range(r), xlab=xlab, ylab=ylab, axes=F, col=col, font.lab=f, lty=lty, cex=cex.mark, ...)
	axis(1, at=r, lwd=3, font=2, cex.axis=cex.axis)
	axis(2, at=c(0, .25, .5, .75, 1), labels=c(0, 25, 50, 75, 100), las=1, lwd=3, font=2, cex.axis=cex.axis)
	if (add.median) names <- paste(names, " (", formatC(tapply(time, group, median, na.rm=T), format="f", digits=1), ")", sep="")
	if (add.surv50) names <- paste(names, " (", formatC(summary(survfit(Surv(time, event) ~ group))$table[,"median"], format="f", digits=1), ")", sep="")
	if (add.legend) legend(legend.pos, names, lty=lty, bty="n", col=col, y.intersp=legend.space, lwd=f, cex=cex.legend)
	#if (add.legend) legend(xrel(.48), yrel(1), names, lty=lty, bty="n", col=col, y.intersp=legend.space, lwd=2, cex=cex.legend)
	if (!is.null(footnote)) mtext(footnote, side=1, line=5, cex=.7)
}





plegend <- function (p, prefix="P") {
	if (p<.000001) paste(prefix, "<0.000001", sep="") else paste(prefix, "=", formatC(p, digits=6, format="f"), sep="")
}


format.mir <- function (m) {
	m <- sub("hsa\\.", "", m)
	m <- gsub("\\.", "-", m)
	m <- sub("mir", "miR", m)
	m <- sub("-*ctrl-*", "/ctrl ", m)
	m <- sub("-$", "\\*", m)
	m
}


pairs <- function (v) {
	do.call(c, sapply(1:(length(v)-1), function (i) lapply((i+1):length(v), function (j) v[c(i, j)])))
}

nemenyi.test <- function(l) {
	t <- oneway_test(values~ind, stack(l), 
			ytrafo = function(data) trafo(data, numeric_trafo = rank),
			xtrafo = function(data) trafo(data, factor_trafo = function(x)
							model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
			teststat = "max", distribution = approximate(B = 10000))
	pvalue(t, method="single-step")
}

nemenyi.pvalue <- function (n, g1, g2) {
	if (paste(g2, "-", g1)%in%rownames(n)) n[paste(g2, "-", g1),] else n[paste(g1, "-", g2),]
}

aov.test <- function (l) {
	aov(values~ind, stack(l))
}

aov.pvalue <- function (a) {
	summary(a)[[1]][1,5]
}

tukey <- function (a, g1, g2) {
	t <- TukeyHSD(a)$ind[,"p adj"]
	if (paste(g2, g1, sep="-")%in%names(t)) t[paste(g2, g1, sep="-")] else t[paste(g1, g2, sep="-")]
}

makeplot <- function (m, l, r, p, normal, legend.pos="topright", ...) {
	name <- paste(format.mir(m), paste(r, collapse=" vs "))
	pointplot(l, normal=T, ylab=format.mir(m), ...)
	legend(legend.pos, plegend(p), bty="n")
	savePlot(file=paste(gsub("/", "", name), ".jpg", sep=""), type="jpeg")
}

read.splitting <- function (s, i, names=c("TRAINING", "VALIDATION")) {
	factor(names[s[,i]])
}



corel.test <- function (u, v) {
	m <- if (min(shapiro.test(u)$p.value, shapiro.test(v)$p.value)<.05) "spearman" else "pearson"
	cor.test(u, v, method="spearman")
}


corelplot <- function (u, v, legend.pos="topleft", col="blue", bold=F, cex.legend=2.2, cex.axis=2.2, ylab="", ylab.adj=0, ...) {
	z <- corel.test(u, v)
	r <- formatC(z$estimat, digits=2, format="f")
	f <- if (bold) 2 else 1
	par(font=f)
	#par(mar=c(5,4+ylab.adj,4,2)+.1)
	par(mar=c(5,6+ylab.adj,4,2)+.1)
	plot(u, v, col=col, xlim=range(pretty(u)), ylim=range(pretty(v)), font.lab=f, axes=F, bty="n", ylab="", ...)
	mtext(ylab, side=2, line=3.5, cex=cex.axis)
	#mtext(ylab, side=2, line=4+ylab.adj, cex=cex.axis)
	#axis(1, lwd=f, font.axis=f, cex.axis=cex.axis)
	#axis(2, lwd=f, las=1, font.axis=f, cex.axis=cex.axis)
	axis(1, lwd=0, lwd.ticks=3, font=2, cex.axis=2.2, cex.lab=2.2)
	axis(2, lwd=0, lwd.ticks=3, las=1, font=2, cex.axis=2.2, cex.lab=2.2, at=pretty(v))
	box(lwd=f)
	legend(legend.pos, legend=paste("r=", r, "  \n", plegend(z$p.value), sep=""), cex=1.4, bty="n")
}

nsplit <- function (v, n, names=1:n) {
	cut(v, breaks=quantile(v, probs=(0:n)/n), include.lowest=T, labels=names)
}

xrel <- function (x) par("usr")[1]+x*diff(par("usr")[1:2])

yrel <- function (y) par("usr")[3]+y*diff(par("usr")[3:4])

myboxplot <- function(l, bold=F, cex.axis=3, footnote=NULL, ylab.adj=0.2, ylab="", group.names=names(l), ...) {
	par(font=2)
	if (!is.null(footnote)) par(oma=c(2,0,0,0))
	r <- pretty(do.call(c, l))
	#par(mar=c(5,4+ylab.adj,4,2)+.1)
	#par(mar=c(10,6,4,2)+.1, font=2, font.lab=2)
	f <- if (bold) 2 else 1
	boxplot(l, axis=F, bty="n", ylab="", ylim=range(r), axes=F, ...)
	#axis(1, lwd=1, font.axis=f, cex.axis=cex.axis, at=seq(l), labels=group.names)
	#axis(2, lwd=f, las=1, font.axis=f, cex.axis=cex.axis, at=r)
	#axis(1, lwd=3, font=2, cex.axis=1.4, cex.lab=1.4)
	#axis(2, lwd=1, lwd.ticks=3, las=1, font=2, cex.axis=1.4, cex.lab=1.4, at=r)
	axis(2, lwd=1, lwd.ticks=3, las=1, font=2, cex.axis=3, cex.lab=3, at=r)#forgrant
	mtext(ylab, side=2, line=3+ylab.adj, cex=cex.axis)
	mtext(group.names, side=1, line=1, at=seq(l), cex=3, font=2)
	mtext(sapply(l, function (y) paste("n=", length(y), sep="")), side=1, line=3.5, at=seq(l), cex=cex.axis, font=2)
	#mtext(sapply(l, function (y) paste("n=", length(which(!is.na(y))), sep="")), side=1, line=3, at=seq(l), cex=cex.axis, font=2)
	box(lwd=3)
	if (!is.null(footnote)) mtext(footnote, side=1, line=5, cex=.7)
}



