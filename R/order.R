
#### Quickselect and Quicksort
## ----------------------------

qselect <- function(x, k = (length(x) + 1L) %/% 2L)
{
	if ( any(k < 1L | k > length(x)) )
		stop("k includes out of bounds subscripts")
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

