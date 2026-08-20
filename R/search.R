
#### Binary search
## ----------------

bsearch <- function(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_)
{
	if ( is.double(query) && is.integer(table) )
		table <- as.double(table)
	if ( is.integer(query) && is.double(table) )
		query <- as.double(query)
	if ( is.unsorted(table) )
		stop("'table' must be sorted")
	relative <- isTRUE(relative)
	ref_side <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_bsearch, query, table, as.double(tolerance),
		relative, ref_side, as.integer(nomatch)) + 1L
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
	if ( ncol(query) != ncol(ref$data) )
		stop("'query' must have the same number of columns as 'ref'")
	if ( anyNA(tolerance) )
		stop("'tolerance' must not contain NAs")
	if ( anyNA(relative) )
		stop("'relative' must not contain NAs")
	tolerance <- as.double(rep_len(tolerance, ncol(ref$data)))
	relative <- ifelse(rep_len(relative, ncol(ref$data)), 1L, 0L)
	.Call(C_do_kdtree_range_search,
		query, ref, tolerance, relative, as.integer(num.threads))
}

