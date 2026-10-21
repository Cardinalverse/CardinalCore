
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

resolve_dx <- function(abs, ppm = numeric(), names = NULL)
{
	k <- max(length(abs), length(ppm), length(names))
	if ( is.list(abs) )
		abs <- unlist(abs)
	if ( is.list(ppm) )
		ppm <- unlist(ppm)
	dx <- rep(abs, length.out=k)
	relative <- setNames(logical(k), names(dx))
	if ( is.null(names(dx)) ) {
		names(dx) <- names
		names(relative) <- names
	} else {
		if ( is.null(names) )
			names <- union(names(dx), names(ppm))
		unexpected <- names(dx) %notin% names
		if ( any(unexpected) ) {
			badnames <- unique(names(dx)[unexpected])
			stop("unexpected name(s) [", paste0(badnames, collapse=", "), "]")
		}
		unresolved <- names %notin% names(dx)
		if ( any(unresolved) ) {
			badnames <- unique(names[unresolved])
			stop("unresolved name(s) [", paste0(badnames, collapse=", "), "]")
		}
		dx <- dx[names]
		relative <- relative[names]
	}
	if ( is.null(names(ppm)) ) {
		i <- seq_along(ppm)
	} else {
		i <- names(ppm)
	}
	dx[i] <- 1e-6 * ppm
	relative[i] <- TRUE
	list(dx=dx, relative=relative, names=names)
}

match_to_self <- function(x, score, tolerance = 0, ppm = numeric(), ...)
{
	x <- as.matrix(x)
	tol <- resolve_dx(tolerance, ppm, colnames(x))
	# match to self reference
	rows <- seq_len(nrow(x))
	hits <- kdsearch(x, x, tolerance=tol$dx,
		relative=tol$relative, relative_to="table")
	hits <- unname(split(hits$index, rep.int(rows, hits$count)))
	# resolve matches based on score
	grp <- vapply(hits, function(i) i[which.max(score[i])], integer(1L))
	# valid hits that match to self-hits
	selfHit <- grp == rows
	validHit <- selfHit[grp]
	# loop until no invalid hits
	while ( any(!validHit) )
	{
		cand <- lapply(hits, intersect, which(selfHit | !validHit))
		grp <- vapply(cand, function(i) i[which.max(score[i])], integer(1L))
		selfHit <- grp == rows
		validHit <- selfHit[grp]
	}
	grp
}

match_to_ref <- function(x, ref, tolerance = 0, ppm = numeric(), ...)
{
	x <- as.matrix(x)
	tol <- resolve_dx(tolerance, ppm, colnames(x))
	if ( is.list(ref) ) {
		# match to sparse reference
		if ( !is.null(names(ref)) && !is.null(tol$names) ) {
			if ( .valid_ref_names(ref, tol$names) )
				ref <- ref[tol$names]
		}
		grp <- msearch(x, ref, tolerance=tol$dx,
			relative=tol$relative, relative_to="table")
		dims <- lengths(ref)
		# convert per-dimension groups to a single linear group
		if ( length(dims) > 1L ) {
			strides <- cumprod(c(1L, dims[-length(dims)]))
			grp <- ((grp - 1L) %*% strides) + 1L
		}
		grp <- as.vector(grp)
	} else {
		ref <- as.matrix(ref)
		if ( !is.null(names(ref)) && !is.null(tol$names) ) {
			if ( .valid_ref_names(ref, tol$names) )
				ref <- ref[tol$names]
		}
		# match to dense reference
		grp <- kdsearch_first(x, ref, tolerance=tol$dx,
			relative=tol$relative, relative_to="table")
	}
	grp
}

.valid_ref_names <- function(ref, names)
{
	if ( !setequal(names(ref), names) )
		stop("invalid 'ref' names [", paste0(names(ref), collapse=", "),
			"]; expected [", paste0(names, collapse=", "), "]")
	TRUE
}

#### Peak pick methods
## -------------------

.peakPick <- function(object, x = "mz", y = "intensity",
	cols = c("snr", "sum", "area", "width", "centroid"), ...)
{
	peaks <- peaks_summary(object[,y,drop=TRUE], object[,x,drop=TRUE], ...)
	cbind(object[peaks$index,,drop=FALSE], do.call(cbind, peaks[cols]))
}

.peakPick_df <- function(object, x = "mz", y = "intensity",
	cols = c("snr", "sum", "area", "width", "centroid"), ...)
{
	peaks <- peaks_summary(object[[y]], object[[x]], ...)
	cbind(object[peaks$index,,drop=FALSE], as.data.frame(peaks[cols]))
}

setMethod("peakPick", "matrix",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakPick(object, x=x, y=y, ...)
})

setMethod("peakPick", "data.frame",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakPick_df(object, x=x, y=y, ...)
})

setMethod("peakPick", "DataFrame",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakPick_df(object, x=x, y=y, ...)
})

#### Peak group methods
## --------------------

.peakGroup <- function(object, x = "mz", y = "intensity", name = "group", ...)
{
	group <- match_to_self(object[,x,drop=FALSE], object[,y,drop=TRUE], ...)
	do.call(cbind, c(list(object), setNames(list(group), name)))
}

.peakGroup_df <- function(object, x = "mz", y = "intensity", name = "group", ...)
{
	group <- match_to_self(object[,x,drop=FALSE], object[[y]], ...)
	do.call(cbind, c(list(object), setNames(list(group), name)))
}

setMethod("peakGroup", "matrix",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakGroup(object, x=x, y=y, ...)
})

setMethod("peakGroup", "data.frame",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakGroup_df(object, x=x, y=y, ...)
})

setMethod("peakGroup", "DataFrame",
	function(object, x = "mz", y = "intensity", ...)
{
	.peakGroup_df(object, x=x, y=y, ...)
})

#### Peak align methods
## --------------------

.peakAlign <- function(object, ref, name = "ref", ...)
{
	if ( is.list(ref) ) {
		x <- names(ref)
	} else {
		x <- colnames(ref)
	}
	if ( is.null(x) )
		stop("'ref' must be a named list or matrix")
	group <- match_to_ref(object[,x,drop=FALSE], ref=ref, ...)
	do.call(cbind, c(list(object), setNames(list(group), name)))
}

setMethod("peakAlign", "matrix",
	function(object, ref, ...)
{
	.peakAlign(object, ref=ref, ...)
})

setMethod("peakAlign", "data.frame",
	function(object, ref, ...)
{
	.peakAlign(object, ref=ref, ...)
})

setMethod("peakAlign", "DataFrame",
	function(object, ref, ...)
{
	.peakAlign(object, ref=ref, ...)
})

#### Peak merge methods
## --------------------

.peakMerge <- function(object, by, count = "n", ...)
{
	if ( by %notin% colnames(object) )
		stop("missing column [", by, "]")
	group <- object[,by,drop=TRUE]
	if ( count %in% colnames(object) ) {
		nobs <- object[,count,drop=TRUE]
	} else {
		nobs <- 1L
	}
	cols <- setdiff(colnames(object), c(by, count))
	stats <- lapply(cols, function(nm) {
		s <- stream_stats(object[,nm,drop=FALSE], nobs=nobs, stat="mean")
		group_stats(s, group=group)
	})
	do.call(cbind, setNames(stats, cols))
}

setMethod("peakMerge", "matrix",
	function(object, by, count = "n", ...)
{
	.peakMerge(object, by=by, count=count, ...)
})


