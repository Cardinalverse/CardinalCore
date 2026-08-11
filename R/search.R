
#### Binary search
## ----------------

bsearch <- function(
	query,
	x,
	tolerance = 0,
	relative = FALSE,
	nearest = FALSE,
	nomatch = NA_integer_)
{
	if ( is.double(query) && is.integer(x) )
		x <- as.double(x)
	if ( is.integer(query) && is.double(x) )
		query <- as.double(query)
	if ( is.unsorted(x) )
		stop("'x' must be sorted")
	.Call(C_do_bsearch, query, x, tolerance,
		isTRUE(relative), isTRUE(nearest), as.integer(nomatch)) + 1L
}

kdtree <- function(data)
{
	if ( inherits(data, "kdtree") )
		return(data)
	if ( is.list(data) )
		data <- do.call(cbind, data)
	data <- as.matrix(data)
	.Call(C_do_kdtree_build, data)
}

kdsearch <- function(
	query,
	data,
	tolerance = 0,
	relative = FALSE,
	num.threads = 1)
{
	if ( !inherits(data, "kdtree") )
		data <- kdtree(data)
	if ( is.null(dim(query)) )
		query <- t(query)
	if ( is.integer(query) && is.double(data$data) )
		storage.mode(query) <- "double"
	if ( is.double(query) && is.integer(data$data) )
		storage.mode(data$data) <- "double"
	if ( is.null(dim(query)) && length(query) != ncol(data$data) )
		stop("query must have the same number of columns as data")
	if ( anyNA(tolerance) )
		stop("tolerance must not contain NAs")
	if ( anyNA(relative) )
		stop("relative must not contain NAs")
	tolerance <- rep_len(as.double(tolerance), ncol(data$data))
	relative <- rep_len(as.logical(relative), ncol(data$data))
	.Call(C_do_kdtree_range_search,
		query, data, tolerance, relative, as.integer(num.threads))
}

