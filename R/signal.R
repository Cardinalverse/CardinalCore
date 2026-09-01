
#### Filter 1D signals
## -------------------

roll <- function(y, k = 5L)
{
	hk <- floor(k / 2)
	lapply(seq_along(y),
		function(i) {
			j <- (i - hk):(i + hk)
			j <- pmin(pmax(j, 1L), length(y))
			y[j]
		})
}

roll_apply <- function(y, k = 5L, FUN = NULL, ..., simplify = TRUE)
{
	y <- lapply(roll(y, k), FUN, ...)
	if ( simplify )
		y <- simplify2array(y)
	y
}

filt1_mean <- function(y, k = 5L)
{
	.Call(C_do_filt1_mean, y, as.integer(k))
}

filt1_conv <- function(y, w = rep(1, 5L))
{
	.Call(C_do_filt1_conv, y, as.double(w))
}


