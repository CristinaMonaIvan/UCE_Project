
pie3D.labels.2 <- function (radialpos, radius = 1, height = 0.1, theta = pi/6, 
		labels, labelcol = par("fg"), labelcex = 1.5, labelrad = NULL, 
		minsep = 0.3) 
{
	if (is.null(labelrad)) labelrad <- rep(1.25, length(labels))
	oldcex <- par("cex")
	nlab <- length(labels)
	par(cex = labelcex, xpd = TRUE)
	for (i in 1:nlab) {
		if (i < nlab) {
			labelsep <- radialpos[i + 1] - radialpos[i]
			if (labelsep < minsep) {
				radialpos[i] <- radialpos[i] + (labelsep - minsep)/2
				radialpos[i + 1] <- radialpos[i + 1] - (labelsep - 
							minsep)/2
			}
		}
		xpos <- labelrad[i] * radius * cos(radialpos[i])
		offset <- (radialpos[i] > pi && radialpos[i] < 2 * pi) * 
				height
		ypos <- labelrad[i] * radius * sin(radialpos[i]) * 2 * theta/pi + 
				sin(radialpos[i]) * height
		text(xpos, ypos, labels[i], col = labelcol, adj = c(0.5, 
						abs(0.5 - sin(radialpos[i])/2)))
	}
	par(cex = oldcex, xpd = FALSE)
}

pie3D <- function (x, edges = NA, radius = 1, height = 0.1, theta = pi/6, 
		start = 0, border = par("fg"), col = NULL, labels = NULL, 
		labelpos = NULL, labelcol = par("fg"), labelcex = 1.5, sector.order = NULL, 
		explode = 0, shade = 0.8, mar = c(4, 4, 4, 4), pty = "s", labelrad=NULL,
		...) 
{
	if (!is.numeric(x) || any(x < 0)) 
		stop("pie3D: x values must be positive numbers")
	if (any(is.na(x))) 
		x <- x[!is.na(x)]
	oldmar <- par("mar")
	par(pty = pty, mar = mar, xpd = TRUE)
	x <- c(0, cumsum(x)/sum(x)) * 2 * pi + start
	nsectors <- length(x) - 1
	if (is.null(col)) 
		col <- rainbow(nsectors)
	else if (length(col) < nsectors) 
		col <- rep(col, length.out = nsectors)
	if (is.null(sector.order)) 
		sector.order <- order(sin((x[2:(nsectors + 1)] + x[1:nsectors])/2), 
				decreasing = TRUE)
	bc <- rep(0, nsectors)
	plot(0, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 
					1), type = "n", axes = FALSE, ...)
	for (i in sector.order) {
		bc[i] <- draw.tilted.sector(radius = radius, height = height, 
				theta = theta, start = x[i], end = x[i + 1], edges = edges, 
				border = border, col = col[i], explode = explode, 
				shade = shade)
	}
	if (!is.null(labels)) {
		if (!is.null(labelpos)) 
			bc <- labelpos
		pie3D.labels.2(bc, height = height, theta = theta, labels = labels, 
				labelcol = labelcol, labelcex = labelcex, labelrad = labelrad)
	}
	par(mar = oldmar, xpd = FALSE, pty = "m")
	invisible(bc)
}