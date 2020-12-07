
split.scale <- function (z, heights, space) {
	x0 <- z[,1]
	y0 <- sapply(1:nrow(z), function (i) if (i==1) 0 else sum(heights[1:(i-1)])+(i-1)*space)
	m <- sapply(1:nrow(z), function (i) heights[i]/(z[i,2]-z[i,1]))
	function (v) {
		if (length(v)==0) return(v)
		i <- sapply(v, function (x) which(sapply(1:nrow(z), function (i) z[i,1]<=x & x<=z[i,2]))[1])
		sapply(seq_along(v), function (j) if (is.na(i[j])) NA else y0[i[j]]+m[i[j]]*(v[j]-x0[i[j]]))
	}
}
