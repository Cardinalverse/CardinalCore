
#### Binary search
## ----------------

bsearch <- function(
	query,
	ref,
	tolerance = 0,
	relative = FALSE,
	nearest = FALSE,
	nomatch = NA_integer_)
{
	if ( is.double(query) && is.integer(ref) )
		ref <- as.double(ref)
	if ( is.integer(query) && is.double(ref) )
		query <- as.double(query)
	if ( is.unsorted(ref) )
		stop("'ref' must be sorted")
	.Call(C_do_bsearch, query, ref, tolerance,
		isTRUE(relative), isTRUE(nearest), as.integer(nomatch)) + 1L
}

kdtree <- function(data)
{
	if ( inherits(data, "kdtree") )
		return(data)
	if ( is.null(dim(data)) ) {
		data <- t(data)
	} else {
		data <- as.matrix(data)
	}
	.Call(C_do_kdtree_build, data)
}

kdsearch <- function(
	query,
	ref,
	tolerance = 0,
	relative = FALSE,
	num.threads = 1)
{
	if ( !inherits(ref, "kdtree") )
		ref <- kdtree(ref)
	if ( is.null(dim(query)) ) {
		query <- t(query)
	} else {
		query <- as.matrix(query)
	}
	if ( is.integer(query) && is.double(ref$data) )
		storage.mode(query) <- "double"
	if ( is.double(query) && is.integer(ref$data) )
		storage.mode(ref$data) <- "double"
	if ( is.null(dim(query)) && length(query) != ncol(ref$data) )
		stop("query must have the same number of columns as data")
	if ( anyNA(tolerance) )
		stop("tolerance must not contain NAs")
	if ( anyNA(relative) )
		stop("relative must not contain NAs")
	tolerance <- as.double(rep_len(tolerance, ncol(ref$data)))
	relative <- ifelse(rep_len(relative, ncol(ref$data)), 1L, 0L)
	.Call(C_do_kdtree_range_search,
		query, ref, tolerance, relative, as.integer(num.threads))
}

