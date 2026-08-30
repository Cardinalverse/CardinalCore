
#### Peak detection
## ----------------

peaks_find <- function(x, k = 5L)
{
	.Call(C_do_peaks_find, x, as.integer(k))
}

peaks_sums <- function(x, k = 5L)
{
	.Call(C_do_peaks_sums, x, as.integer(k))
}


peaks_prominences <- function(x, k = 5L, wlen = length(x))
{
	.Call(C_do_peaks_prominences, x, as.integer(k), as.integer(wlen))
}


