
#### Compute differences
## ---------------------

qdiff <- function(x, ref, relative = FALSE)
{
	if ( missing(ref) ) {
		ref <- x[-length(x)]
		x <- x[-1L]
		len <- length(x)
	} else {
		len <- max(length(x), length(ref))
	}
	if ( is.integer(x) && is.double(ref) )
		x <- as.double(x)
	if ( is.double(x) && is.integer(ref) )
		ref <- as.double(ref)
	x <- rep_len(x, len)
	ref <- rep_len(ref, len)
	.Call(C_do_qdiff, x, ref, isTRUE(relative))
}

#### Quickselect and Quicksort
## ----------------------------

qselect <- function(x, k = (length(x) + 1L) %/% 2L)
{
	if ( any(k < 1L | k > length(x)) )
		stop("k includes out of bounds indices")
	.Call(C_do_qselect, x, as.integer(k - 1L))
}

qorder <- function(x)
{
	.Call(C_do_qorder, x) + 1L
}

#### Median and MAD
## -----------------

qmedian <- function(x)
{
	.Call(C_do_qmedian, x)
}

qmad <- function(x, center = qmedian(x), constant = 1.4826)
{
	.Call(C_do_qmad, x, center, constant)
}

#### Binary search
## ----------------

bsearch <- function(x, data, tolerance = 0,
	relative = FALSE, nearest = FALSE, nomatch = NA_integer_)
{
	if ( is.integer(x) && is.double(data) )
		x <- as.double(x)
	if ( is.double(x) && is.integer(data) )
		data <- as.double(data)
	if ( is.unsorted(data) )
		stop("'data' must be sorted")
	.Call(C_do_bsearch, x, data, tolerance,
		isTRUE(relative), isTRUE(nearest), as.integer(nomatch)) + 1L
}
