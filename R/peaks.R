
#### Peak detection
## ----------------

peaks_find <- function(x, k = 5L)
{
	.Call(C_do_peaks_find, x, as.integer(k))
}

peaks_find_limits <- function(x, k = 5L)
{
	.Call(C_do_peaks_find_limits, x, as.integer(k))
}



