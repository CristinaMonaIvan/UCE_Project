alluvial <- function (..., freq, col = "gray", border = 0, layer, hide = FALSE, 
		alpha = 0.5, gap.width = 0.05, xw = 0.1, cw = 0.1, blocks = TRUE, 
		ordering = NULL, axis_labels = NULL, cex = par("cex"), cex.axis = par("cex.axis"),
		rect.col=NULL) 
{
	p <- data.frame(..., freq = freq, col, alpha, border, hide, 
			stringsAsFactors = FALSE)
	np <- ncol(p) - 5
	if (!is.null(ordering)) {
		stopifnot(is.list(ordering))
		if (length(ordering) != np) 
			stop("'ordering' argument should have ", np, " components, has ", 
					length(ordering))
	}
	n <- nrow(p)
	if (missing(layer)) {
		layer <- 1:n
	}
	p$layer <- layer
	d <- p[, 1:np, drop = FALSE]
	p <- p[, -c(1:np), drop = FALSE]
	p$freq <- with(p, freq/sum(freq))
	col <- col2rgb(p$col, alpha = TRUE)
	if (!identical(alpha, FALSE)) {
		col["alpha", ] <- p$alpha * 256
	}
	p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), 
								maxColorValue = 256)))
	isch <- sapply(d, is.character)
	d[isch] <- lapply(d[isch], as.factor)
	if (length(blocks) == 1) {
		blocks <- if (!is.na(as.logical(blocks))) {
					rep(blocks, np)
				}
				else if (blocks == "bookends") {
					c(TRUE, rep(FALSE, np - 2), TRUE)
				}
	}
	if (is.null(axis_labels)) {
		axis_labels <- names(d)
	}
	else {
		if (length(axis_labels) != ncol(d)) 
			stop("`axis_labels` should have length ", names(d), 
					", has ", length(axis_labels))
	}
	getp <- function(i, d, f, w = gap.width) {
		#a <- c(i, (1:ncol(d))[-i])
		a <- i
		if (is.null(ordering[[i]])) {
			o <- do.call(order, d[a])
		}
		else {
			#d2 <- d
			#d2[1] <- ordering[[i]]
			o <- ordering[[i]]
			#o <- do.call(order, d2[a])
		}
		x <- c(0, cumsum(f[o])) * (1 - w)
		x <- cbind(x[-length(x)], x[-1])
		gap <- cumsum(c(0L, diff(as.numeric(d[o, i])) != 0))
		mx <- max(gap)
		if (mx == 0) 
			mx <- 1
		gap <- gap/mx * w
		(x + gap)[order(o), ]
	}
	dd <- lapply(seq_along(d), getp, d = d, f = p$freq)
	rval <- list(endpoints = dd)
	op <- par(mar = c(2, 1, 1, 1))
	plot(NULL, type = "n", xlim = c(1 - cw[1], np + cw[np]), ylim = c(0, 
					1), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", 
			ylab = "", frame = FALSE)
	ind <- which(!p$hide)[rev(order(p[!p$hide, ]$layer))]
	for (i in ind) {
		for (j in 1:(np - 1)) {
			xspline(c(j, j, j + xw, j + 1 - xw, j + 1, j + 1, 
									j + 1 - xw, j + xw, j) + rep(c(cw[j], -cw[j+1], cw[j]), 
									c(3, 4, 2)), c(dd[[j]][i, c(1, 2, 2)], rev(dd[[j + 
													1]][i, c(1, 1, 2, 2)]), dd[[j]][i, c(1, 1)]), 
					shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0), open = FALSE, 
					col = p$col[i], border = p$border[i])
		}
	}
	for (j in seq_along(dd)) {
		ax <- lapply(split(dd[[j]], d[, j]), range)
		if (blocks[j]) {
			for (k in seq_along(ax)) {
				rect(j - cw[j], ax[[k]][1], j + cw[j], ax[[k]][2], col=rect.col[[j]][names(ax)[k]])
			}
		}
		else {
			for (i in ind) {
				x <- j + c(-1, 1) * cw
				y <- t(dd[[j]][c(i, i), ])
				w <- xw * (x[2] - x[1])
				xspline(x = c(x[1], x[1], x[1] + w, x[2] - w, 
								x[2], x[2], x[2] - w, x[1] + w, x[1]), y = c(y[c(1, 
												2, 2), 1], y[c(2, 2, 1, 1), 2], y[c(1, 1), 
										1]), shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0), 
						open = FALSE, col = p$col[i], border = p$border[i])
			}
		}
		for (k in seq_along(ax)) {
			text(j, mean(ax[[k]]), labels = names(ax)[k], cex = cex[j])
		}
	}
	#axis(1, at = c(rbind(cw, -cw)) + rep(seq_along(d), 
	#				each = 2), line = 0.5, col = "white", col.ticks = "black", 
	#		labels = FALSE)
	#axis(1, at = seq_along(d), tick = FALSE, labels = axis_labels, 
	#		cex.axis = cex.axis)
	par(op)
	invisible(rval)
}
