
#### Peak summarization
## --------------------

peaks_find <- function(y, k = 5L)
{
	k <- as.integer(min(max(3L, k), length(y)))
	.Call(C_do_peaks_find, y, as.integer(k))
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
	noise = c("DiffMAD", "SmoothSD", "SmoothMAD"), wlen = 0L, fmax = 0.5)
{
	k <- as.integer(min(max(3L, k), length(y)))
	fmax <- as.double(min(max(fmax, 0), 1))
	wlen <- as.integer(min(max(0L, wlen), length(y)))
	noise <- c("DiffMAD"=0L, "SmoothSD"=1L, "SmoothMAD"=2L)[match.arg(noise)]
	.Call(C_do_peaks_summary, y, x, k, noise, wlen, fmax)
}

peaks_agg <- function(y, x = seq_along(y), xout,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0, ppm = numeric())
{
	if ( is.unsorted(x) )
	{
		i <- order(x)
		y <- y[i]
		x <- x[i]
	}
	relative <- !is.na(ppm)
	if ( !is.na(ppm) )
		tolerance <- 1e-6 * ppm
	bsearch_agg(xout, x, y, stat=stat, tolerance=tolerance, relative=relative)
}

#### Peak alignment
## --------------------

resolve_dx <- function(abs, ppm, names)
{
	dx <- rep(abs, length.out=length(names))
	relative <- logical(length(names))
	if ( is.null(names(dx)) ) {
		names(dx) <- names
		names(relative) <- names
	}
	if ( is.null(names(ppm)) ) {
		i <- seq_along(ppm)
	} else {
		i <- names(ppm)
	}
	dx[i] <- 1e-6 * ppm
	relative[i] <- TRUE
	nomatch <- names(dx) %notin% names
	if ( any(nomatch) ) {
		badnames <- names(dx)[nomatch]
		stop("unexpected name(s) [", paste0(badnames, collapse=", "), "]")
	}
	list(dx=dx, relative=relative)
}

peakslist_match <- function(peakslist, ref, stat = "mean",
	tolerance = 0, ppm = numeric())
{
	nrs <- vapply(peakslist, NROW, numeric(1L))
	
}

