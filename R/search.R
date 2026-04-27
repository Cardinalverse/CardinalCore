
sdiff <- function(x, ref, units = c("absolute", "relative", "ppm"))
{
	units <- match.arg(units)
	if ( missing(ref) ) {
		ref <- x[-1L]
		x <- x[-length(x)]
	}
	if ( is.integer(x) && is.double(ref) )
		x <- as.double(x)
	if ( is.double(x) && is.integer(ref) )
		y <- as.double(ref)
	.Call(C_compute_sdiff, x, ref, units != "absolute")
}
