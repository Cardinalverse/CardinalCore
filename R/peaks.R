
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
	if ( is.unsorted(x) )
		stop("'x' must be sorted")
	k <- as.integer(min(max(3L, k), length(y)))
	fmax <- as.double(min(max(fmax, 0), 1))
	.Call(C_do_peaks_widths, y, x, k, fmax)
}

peaks_areas <- function(y, x = seq_along(y), k = 5L)
{
	if ( is.unsorted(x) )
		stop("'x' must be sorted")
	k <- as.integer(min(max(3L, k), length(y)))
	.Call(C_do_peaks_areas, y, x, k)
}

peaks_summary <- function(y, x = seq_along(y), k = 5L, 
	noise = c("DiffMAD", "SmoothSD", "SmoothMAD"),
	wlen = 0L, fmax = 0.5, ...)
{
	if ( is.unsorted(x) )
		stop("'x' must be sorted")
	k <- as.integer(min(max(3L, k), length(y)))
	fmax <- as.double(min(max(fmax, 0), 1))
	wlen <- as.integer(min(max(0L, wlen), length(y)))
	noise <- c("DiffMAD"=0L, "SmoothSD"=1L, "SmoothMAD"=2L)[match.arg(noise)]
	.Call(C_do_peaks_summary, y, x, k, noise, wlen, fmax)
}

#### Peak grouping
## ----------------

resolve_dx <- function(abs, ppm, names)
{
	k <- max(length(abs), length(ppm), length(names))
	dx <- rep(abs, length.out=k)
	relative <- logical(k)
	if ( is.null(names(dx)) ) {
		names(dx) <- names
		names(relative) <- names
	}
	if ( is.null(names(ppm)) ) {
		i <- seq_along(ppm)
	} else {
		i <- names(ppm)
	}
	if ( !is.null(names(dx)) && !is.null(names) ) {
		dx <- dx[names]
		relative <- relative[names]
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

group_by_ref <- function(x, ref, tolerance = 0, ppm = numeric(), ...)
{
	if ( is.numeric(ref) && !is.array(ref) )
		ref <- list(ref)
	if ( is.list(ref) ) {
		tol <- resolve_dx(tolerance, ppm, names(ref))
		if ( !identical(colnames(x), names(ref)) ) {
			if ( is.null(colnames(x)) || is.null(names(ref)) ) {
				x <- x[,seq_along(ref),drop=FALSE]
			} else {
				x <- x[,names(ref),drop=FALSE]
			}
		}
		grp <- msearch(x, ref, tolerance=tol$dx,
			relative=tol$relative, relative_to="table")
		dims <- lengths(ref)
		if ( length(dims) > 1L ) {
			strides <- cumprod(c(1L, dims[-length(dims)]))
			grp <- ((grp - 1L) %*% strides) + 1L
		}
		grp <- as.vector(grp)
	} else if ( is.matrix(ref) ) {
		tol <- resolve_dx(tolerance, ppm, colnames(ref))
		if ( !identical(colnames(x), colnames(ref)) ) {
			if ( is.null(colnames(x)) || is.null(colnames(ref)) ) {
				x <- x[,seq_len(ncol(ref)),drop=FALSE]
			} else {
				x <- x[,colnames(ref),drop=FALSE]
			}
		}
		grp <- kdsearch_first(x, ref, tolerance=tol$dx,
			relative=tol$relative, relative_to="table")
	} else {
		stop("'ref' must be a numeric vector, list or matrix")
	}
	grp
}

#### Peak pick methods
## -------------------

.peakPicked <- function(object, peaks,
	cols = c("snr", "sum", "area", "width", "centroid"))
{
	object <- object[peaks$index,,drop=FALSE]
	do.call(cbind, c(list(object), peaks[cols]))
}

setMethod("peakPick", "matrix",
	function(object, x = "mz", y = "intensity", ...)
{
	peaks <- peaks_summary(object[,y], object[,x], ...)
	as.matrix(.peakPicked(object, peaks))
})

setMethod("peakPick", "data.frame",
	function(object, x = "mz", y = "intensity", ...)
{
	peaks <- peaks_summary(object[[y]], object[[x]], ...)
	as.data.frame(.peakPicked(object, peaks))
})

setMethod("peakPick", "DataFrame",
	function(object, x = "mz", y = "intensity", ...)
{
	peaks <- peaks_summary(object[[y]], object[[x]], ...)
	DataFrame(.peakPicked(object, peaks))
})

#### Peak align methods
## --------------------

.peakAligned <- function(object, group, name = "ref")
{
	grouplist <- setNames(list(group), name)
	do.call(cbind, c(list(object), grouplist))
}

setMethod("peakAlign", "matrix",
	function(object, ref, ...)
{
	group <- group_by_ref(object, ref, ...)
	as.matrix(.peakAligned(object, group))
})

setMethod("peakAlign", "data.frame",
	function(object, ref, ...)
{
	group <- group_by_ref(object, ref, ...)
	as.data.frame(.peakAligned(object, group))
})

setMethod("peakAlign", "DataFrame",
	function(object, ref, ...)
{
	group <- group_by_ref(object, ref, ...)
	DataFrame(.peakAligned(object, group))
})


