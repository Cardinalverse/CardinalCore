
#### Peak detection
## ----------------

peaks_find <- function(y, k = 5L)
{
	k <- as.integer(min(max(3L, k), length(y)))
	.Call(C_do_peaks_find, y, as.integer(k))
}

peaks_snrs <- function(y, k = 5L,
	noise = c("Diff", "SmoothSD", "SmoothMAD"), wlen = 0L)
{
	k <- as.integer(min(max(3L, k), length(y)))
	wlen <- as.integer(min(max(0L, wlen), length(y)))
	noise <- c("Diff"=0L, "SmoothSD"=1L, "SmoothMAD"=2L)[match.arg(noise)]
	.Call(C_do_peaks_snrs, y, k, noise, wlen)
}

peaks_prominences <- function(y, k = 5L, wlen = 0L)
{
	k <- as.integer(min(max(3L, k), length(y)))
	wlen <- as.integer(min(max(0L, wlen), length(y)))
	.Call(C_do_peaks_prominences, y, k, wlen)
}

peaks_widths <- function(y, x = seq_along(y), k = 5L, fmax = 0.5)
{
	k <- as.integer(min(max(3L, k), length(y)))
	fmax <- as.double(min(max(fmax, 0), 1))
	.Call(C_do_peaks_widths, y, x, k, fmax)
}

peaks_areas <- function(y, x = seq_along(y), k = 5L)
{
	k <- as.integer(min(max(3L, k), length(y)))
	.Call(C_do_peaks_areas, y, x, k)
}

peaks_summary <- function(y, x = seq_along(y), k = 5L, 
	noise = c("Diff", "SmoothSD", "SmoothMAD"), wlen = 0L, fmax = 0.5)
{
	k <- as.integer(min(max(3L, k), length(y)))
	fmax <- as.double(min(max(fmax, 0), 1))
	wlen <- as.integer(min(max(0L, wlen), length(y)))
	noise <- c("Diff"=0L, "SmoothSD"=1L, "SmoothMAD"=2L)[match.arg(noise)]
	.Call(C_do_peaks_summary, y, x, k, noise, wlen, fmax)
}

