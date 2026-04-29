
#### Compute differences
## ---------------------

qdiff <- function(x, ref, units = c("absolute", "relative", "ppm"))
{
	units <- match.arg(units)
	if ( missing(ref) ) {
		ref <- x[-length(x)]
		x <- x[-1L]
	}
	if ( is.integer(x) && is.double(ref) )
		x <- as.double(x)
	if ( is.double(x) && is.integer(ref) )
		ref <- as.double(ref)
	len <- max(length(x), length(ref))
	x <- rep_len(x, len)
	ref <- rep_len(ref, len)
	q <- switch(units, ppm=1e6, 1)
	q * .Call(C_do_qdiff, x, ref, units != "absolute")
}

#### Quickselect and Quicksort
## ----------------------------

qselect <- function(x, k = (length(x) + 1L) %/% 2L)
{
	if ( any(k < 1L | k > length(x)) )
		stop("k includes out of bounds indices")
	.Call(C_do_qselect, x, as.integer(k - 1L))
}

qsort <- function(x, decreasing = FALSE, index.return = FALSE)
{
	ix <- .Call(C_do_qsort, x) + 1L
	if ( decreasing )
		ix <- rev(ix)
	if ( index.return ) {
		list(x=x[ix], ix=ix)
	} else {
		x[ix]
	}
}

#### Median and MAD
## -----------------

qmedian <- function(x, na.rm = FALSE)
{
	if ( !na.rm && anyNA(x) ) {
		NA_real_
	} else {
		.Call(C_do_qmedian, x)
	}
}

qmad <- function(x, center = qmedian(x), constant = 1.4826, na.rm = FALSE)
{
	if ( !na.rm && anyNA(x) ) {
		NA_real_
	} else {
		.Call(C_do_qmad, x, center, constant)
	}
}

