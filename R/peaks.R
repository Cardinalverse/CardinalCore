
#### Peak detection
## ----------------

peaks_find <- function(y, k = 5L)
{
	.Call(C_do_peaks_find, y, as.integer(k))
}

peaks_snrs <- function(y, k = 5L, method = "Diff", wlen = 0)
{
	method <- c("Diff"=0L, "SmoothSD"=1L, "SmoothMAD"=2L)[[method]]
	.Call(C_do_peaks_snrs, y, as.integer(k), method, as.integer(wlen))
}

peaks_sums <- function(y, k = 5L)
{
	.Call(C_do_peaks_sums, y, as.integer(k))
}

peaks_prominences <- function(y, k = 5L, wlen = 0)
{
	.Call(C_do_peaks_prominences, y, as.integer(k), as.integer(wlen))
}

peaks_areas <- function(y, x = seq_along(y), k = 5L)
{
	if ( is.double(y) && is.integer(x) )
		x <- as.double(x)
	if ( is.integer(y) && is.double(x) )
		y <- as.double(y)
	.Call(C_do_peaks_areas, y, x, as.integer(k))
}

peaks_centroids <- function(y, x = seq_along(y), k = 5L)
{
	if ( is.double(y) && is.integer(x) )
		x <- as.double(x)
	if ( is.integer(y) && is.double(x) )
		y <- as.double(y)
	.Call(C_do_peaks_centroids, y, x, as.integer(k))
}

peaks_widths <- function(y, x = seq_along(y), k = 5L, fmax = 0.5)
{
	if ( is.double(y) && is.integer(x) )
		x <- as.double(x)
	if ( is.integer(y) && is.double(x) )
		y <- as.double(y)
	.Call(C_do_peaks_widths, y, x, as.integer(k), as.double(fmax))
}

