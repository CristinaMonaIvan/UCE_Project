mosaicplot <- function (x, main = deparse(substitute(x)), sub = NULL, xlab = NULL, 
		ylab = NULL, sort = NULL, off = NULL, dir = NULL, color = NULL, 
		shade = FALSE, margin = NULL, cex.axis = 0.66, las = par("las"), 
		border = NULL, type = c("pearson", "deviance", "FT"), 
		rect.label=NULL, rect.cex=1, rect.font=1, ...) 
{
	coord <- function (X, x1, y1, x2, y2, off, maxdim) {
		p <- ncol(X)-2
		xdim <- maxdim[1L]
		XP <- rep.int(0, xdim)
		for (i in seq_len(xdim)) XP[i] <- sum(X[X[, 1L] == i, p])/sum(X[, p])
		white <- off[1L] * (x2 - x1)/max(1, xdim - 1)
		x.l <- x1
		x.r <- x1 + (1 - off[1L]) * XP[1L] * (x2 - x1)
		if (xdim > 1L) 
			for (i in 2:xdim) {
				x.l <- c(x.l, x.r[i - 1L] + white)
				x.r <- c(x.r, x.r[i - 1L] + white + (1 - off[1L]) * XP[i] * (x2 - x1))
			}
		Y <- apply(x, 1, function (y) c(y2, sapply(seq(y), function (i) y2+(y1-y2)*sum(y[1:i])/sum(y))))
		y <- apply(Y, 2, function (y) sapply(1:(length(y)-1), function (i) (y[i]+y[i+1])/2))
		list(x=(x.l+x.r)/2, y=y)
	}
	mosaic.cell <- function(X, x1, y1, x2, y2, srt.x, srt.y, 
			adj.x, adj.y, off, dir, color, lablevx, lablevy, maxdim, 
			currlev, label) {
		p <- ncol(X) - 2
		if (dir[1L] == "v") {
			xdim <- maxdim[1L]
			XP <- rep.int(0, xdim)
			for (i in seq_len(xdim)) XP[i] <- sum(X[X[, 1L] == 
										i, p])/sum(X[, p])
			if (anyNA(XP)) 
				stop("missing values in contingency table")
			white <- off[1L] * (x2 - x1)/max(1, xdim - 1)
			x.l <- x1
			x.r <- x1 + (1 - off[1L]) * XP[1L] * (x2 - x1)
			if (xdim > 1L) 
				for (i in 2:xdim) {
					x.l <- c(x.l, x.r[i - 1L] + white)
					x.r <- c(x.r, x.r[i - 1L] + white + (1 - off[1L]) * 
									XP[i] * (x2 - x1))
				}
			if (lablevx > 0L) {
				this.lab <- if (is.null(label[[1L]][1L])) {
							rep.int("", length(currlev))
						}
						else label[[1L]]
				text(x = x.l + (x.r - x.l)/2, y = 1000 - 35 * 
								cex.axis/0.66 + 22 * cex.axis/0.65 * (lablevx - 
									1), srt = srt.x, adj = adj.x, cex = cex.axis, 
						this.lab, xpd = NA)
				
			}
			if (p > 2L) {
				for (i in seq_len(xdim)) {
					if (XP[i] > 0) {
						Recall(X[X[, 1L] == i, 2L:(p + 2L), drop = FALSE], 
								x.l[i], y1, x.r[i], y2, srt.x, srt.y, adj.x, 
								adj.y, off[-1L], dir[-1L], color, lablevx - 
										1, (i == 1L) * lablevy, maxdim[-1L], 
								currlev + 1, label[2:p])
					}
					else {
						segments(rep.int(x.l[i], 3L), y1 + (y2 - 
											y1) * c(0, 2, 4)/5, rep.int(x.l[i], 3L), 
								y1 + (y2 - y1) * c(1, 3, 5)/5)
					}
				}
			}
			else {
				for (i in seq_len(xdim)) {
					if (XP[i] > 0) {
						polygon(c(x.l[i], x.r[i], x.r[i], x.l[i]), 
								c(y1, y1, y2, y2), lty = if (extended) 
											X[i, p + 1L]
										else 1L, col = color[if (extended) 
													X[i, p + 2L]
												else i], border = border)
					}
					else {
						segments(rep.int(x.l[i], 3L), y1 + (y2 - 
											y1) * c(0, 2, 4)/5, rep.int(x.l[i], 3L), 
								y1 + (y2 - y1) * c(1, 3, 5)/5)
					}
				}
			}
		}
		else {
			ydim <- maxdim[1L]
			YP <- rep.int(0, ydim)
			for (j in seq_len(ydim)) {
				YP[j] <- sum(X[X[, 1L] == j, p])/sum(X[, p])
			}
			white <- off[1L] * (y2 - y1)/(max(1, ydim - 1))
			y.b <- y2 - (1 - off[1L]) * YP[1L] * (y2 - y1)
			y.t <- y2
			if (ydim > 1L) {
				for (j in 2:ydim) {
					y.b <- c(y.b, y.b[j - 1] - white - (1 - off[1L]) * 
									YP[j] * (y2 - y1))
					y.t <- c(y.t, y.b[j - 1] - white)
				}
			}
			if (lablevy > 0L) {
				this.lab <- if (is.null(label[[1L]][1L])) {
							rep.int("", length(currlev))
						}
						else label[[1L]]
				text(x = 35 * cex.axis/0.66 - 20 * cex.axis/0.66 * 
								(lablevy - 1), y = y.b + (y.t - y.b)/2, srt = srt.y, 
						adj = adj.y, cex = cex.axis, this.lab, xpd = NA)
			}
			if (p > 2L) {
				for (j in seq_len(ydim)) {
					if (YP[j] > 0) {
						Recall(X[X[, 1L] == j, 2:(p + 2), drop = FALSE], 
								x1, y.b[j], x2, y.t[j], srt.x, srt.y, adj.x, 
								adj.y, off[-1L], dir[-1L], color, (j == 
											1L) * lablevx, lablevy - 1, maxdim[-1L], 
								currlev + 1, label[2:p])
					}
					else {
						segments(x1 + (x2 - x1) * c(0, 2, 4)/5, rep.int(y.b[j], 
										3L), x1 + (x2 - x1) * c(1, 3, 5)/5, rep.int(y.b[j], 
										3L))
					}
				}
			}
			else {
				for (j in seq_len(ydim)) {
					if (YP[j] > 0) {
						polygon(c(x1, x2, x2, x1), c(y.b[j], y.b[j], 
										y.t[j], y.t[j]), lty = if (extended) 
											X[j, p + 1]
										else 1, col = color[if (extended) 
													X[j, p + 2]
												else j], border = border)
					}
					else {
						segments(x1 + (x2 - x1) * c(0, 2, 4)/5, rep.int(y.b[j], 
										3L), x1 + (x2 - x1) * c(1, 3, 5)/5, rep.int(y.b[j], 
										3L))
					}
				}
			}
		}
		invisible()
	}
	srt.x <- if (las > 1) 
				90
			else 0
	srt.y <- if (las == 0 || las == 3) 
				90
			else 0
	if (is.null(dim(x))) 
		x <- as.array(x)
	else if (is.data.frame(x)) 
		x <- data.matrix(x)
	dimd <- length(dx <- dim(x))
	if (dimd == 0L || any(dx == 0L)) 
		stop("'x' must not have 0 dimensionality")
	chkDots(...)
	Ind <- 1L:dx[1L]
	if (dimd > 1L) {
		Ind <- rep.int(Ind, prod(dx[2:dimd]))
		for (i in 2:dimd) {
			Ind <- cbind(Ind, c(matrix(1L:dx[i], byrow = TRUE, 
									nrow = prod(dx[1L:(i - 1)]), ncol = prod(dx[i:dimd]))))
		}
	}
	Ind <- cbind(Ind, c(x))
	if (is.logical(shade) && !shade) {
		extended <- FALSE
		Ind <- cbind(Ind, NA, NA)
	}
	else {
		if (is.logical(shade)) 
			shade <- c(2, 4)
		else if (any(shade <= 0) || length(shade) > 5) 
			stop("invalid 'shade' specification")
		extended <- TRUE
		shade <- sort(shade)
		breaks <- c(-Inf, -rev(shade), 0, shade, Inf)
		color <- c(hsv(0, s = seq.int(1, to = 0, length.out = length(shade) + 
										1)), hsv(4/6, s = seq.int(0, to = 1, length.out = length(shade) + 
										1)))
		if (is.null(margin)) 
			margin <- as.list(1L:dimd)
		E <- stats::loglin(x, margin, fit = TRUE, print = FALSE)$fit
		type <- match.arg(type)
		residuals <- switch(type, pearson = (x - E)/sqrt(E), 
				deviance = {
					tmp <- 2 * (x * log(ifelse(x == 0, 1, x/E)) - 
								(x - E))
					tmp <- sqrt(pmax(tmp, 0))
					ifelse(x > E, tmp, -tmp)
				}, FT = sqrt(x) + sqrt(x + 1) - sqrt(4 * E + 1))
		Ind <- cbind(Ind, c(1 + (residuals < 0)), as.numeric(cut(residuals, 
								breaks)))
	}
	label <- dimnames(x)
	if (is.null(off)) 
		off <- if (dimd == 2) 
					2 * (dx - 1)
				else rep.int(10, dimd)
	if (length(off) != dimd) 
		off <- rep_len(off, dimd)
	if (any(off > 50)) 
		off <- off * 50/max(off)
	if (is.null(dir) || length(dir) != dimd) {
		dir <- rep_len(c("v", "h"), dimd)
	}
	if (!is.null(sort)) {
		if (length(sort) != dimd) 
			stop("length of 'sort' does not conform to 'dim(x)'")
		Ind[, seq_len(dimd)] <- Ind[, sort]
		off <- off[sort]
		dir <- dir[sort]
		label <- label[sort]
	}
	nam.dn <- names(label)
	if (is.null(xlab) && any(dir == "v")) 
		xlab <- nam.dn[min(which(dir == "v"))]
	if (is.null(ylab) && any(dir == "h")) 
		ylab <- nam.dn[min(which(dir == "h"))]
	ncolors <- length(tabulate(Ind[, dimd]))
	if (!extended && ((is.null(color) || length(color) != ncolors))) {
		color <- if (is.logical(color)) 
					if (color[1L]) 
						gray.colors(ncolors)
					else rep.int(0, ncolors)
				else if (is.null(color)) 
					rep.int("grey", ncolors)
				else rep_len(color, ncolors)
	}
	dev.hold()
	on.exit(dev.flush())
	plot.new()
	if (!extended) {
		opar <- par(usr = c(1, 1000, 1, 1000), mgp = c(1, 1, 
						0))
		on.exit(par(opar), add = TRUE)
	}
	else {
		pin <- par("pin")
		rtxt <- "Standardized\nResiduals:"
		rtxtCex <- min(1, pin[1L]/(strheight(rtxt, units = "inches") * 
							12), pin[2L]/(strwidth(rtxt, units = "inches")/4))
		rtxtWidth <- 0.1
		opar <- par(usr = c(1, 1000 * (1.1 + rtxtWidth), 1, 1000), 
				mgp = c(1, 1, 0))
		on.exit(par(opar), add = TRUE)
		rtxtHeight <- strwidth(rtxt, units = "i", cex = rtxtCex)/pin[2L]
		text(1000 * (1.05 + 0.5 * rtxtWidth), 0, labels = rtxt, 
				adj = c(0, 0.25), srt = 90, cex = rtxtCex)
		len <- length(shade) + 1
		bh <- 0.95 * (0.95 - rtxtHeight)/(2 * len)
		x.l <- 1000 * 1.05
		x.r <- 1000 * (1.05 + 0.7 * rtxtWidth)
		y.t <- 1000 * rev(seq.int(from = 0.95, by = -bh, length.out = 2 * 
								len))
		y.b <- y.t - 1000 * 0.8 * bh
		ltype <- c(rep.int(2, len), rep.int(1, len))
		for (i in 1:(2 * len)) {
			polygon(c(x.l, x.r, x.r, x.l), c(y.b[i], y.b[i], 
							y.t[i], y.t[i]), col = color[i], lty = ltype[i], 
					border = border)
		}
		brks <- round(breaks, 2)
		y.m <- y.b + 1000 * 0.4 * bh
		text(1000 * (1.05 + rtxtWidth), y.m, c(paste0("<", brks[2L]), 
						paste(brks[2:(2 * len - 1)], brks[3:(2 * len)], sep = ":"), 
						paste0(">", brks[2 * len])), srt = 90, cex = cex.axis, 
				xpd = NA)
	}
	if (!is.null(main) || !is.null(xlab) || !is.null(ylab) || 
			!is.null(sub)) 
		title(main, sub = sub, xlab = xlab, ylab = ylab)
	adj.x <- adj.y <- 0.5
	x1 <- 30 + 20 * cex.axis/0.66
	y1 <- 5
	x2 <- 950
	y2 <- 1000 - x1
	maxlen.xlabel <- maxlen.ylabel <- 35 * cex.axis/0.66
	if (srt.x == 90) {
		maxlen.xlabel <- max(strwidth(label[[dimd + 1L - match("v", 
										rev(dir))]], cex = cex.axis))
		adj.x <- 1
		y2 <- y2 - maxlen.xlabel
	}
	if (srt.y == 0) {
		maxlen.ylabel <- max(strwidth(label[[match("h", dir)]], 
						cex = cex.axis))
		adj.y <- 0
		x1 <- x1 + maxlen.ylabel
	}
	mosaic.cell(Ind, x1 = x1, y1 = y1, x2 = x2, y2 = y2, srt.x = srt.x, 
			srt.y = srt.y, adj.x = adj.x, adj.y = adj.y, off = off/100, 
			dir = dir, color = color, lablevx = 2, lablevy = 2, maxdim = apply(as.matrix(Ind[, 
									1L:dimd]), 2L, max), currlev = 1, label = label)
	z <- coord(Ind, x1, y1, x2, y2, off/100, maxdim=apply(as.matrix(Ind[,1L:dimd]), 2L, max))
	if (!is.null(rect.label)) {
		z <- coord(Ind, x1, y1, x2, y2, off/100, maxdim=apply(as.matrix(Ind[,1L:dimd]), 2L, max))
		for (i in 1:nrow(x)) {
			for (j in 1:ncol(x))
			{
				if (x[i,j]>0) text(z$x[i], z$y[j,i], rect.label[i,j], cex=rect.cex, font=rect.font)
			}
		}
	}
}