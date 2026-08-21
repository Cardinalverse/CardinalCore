
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
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_bsearch, query, table, as.double(tolerance),
		relative, referent, as.integer(nomatch)) + 1L
}

kdtree <- function(table)
{
	if ( inherits(table, "kdtree") )
		return(table)
	if ( is.null(dim(table)) ) {
		table <- t(table)
	} else {
		table <- as.matrix(table)
	}
	.Call(C_do_kdtree_build, table)
}

kdsearch <- function(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	num.threads = 1)
{
	if ( !inherits(table, "kdtree") )
		table <- kdtree(table)
	if ( is.null(dim(query)) ) {
		query <- t(query)
	} else {
		query <- as.matrix(query)
	}
	if ( is.integer(query) && is.double(table$table) )
		storage.mode(query) <- "double"
	if ( is.double(query) && is.integer(table$table) )
		storage.mode(table$table) <- "double"
	if ( ncol(query) != ncol(table$table) )
		stop("'query' must have the same number of columns as 'table'")
	if ( anyNA(tolerance) )
		stop("'tolerance' must not contain NAs")
	if ( anyNA(relative) )
		stop("'relative' must not contain NAs")
	tolerance <- as.double(rep_len(tolerance, ncol(table$table)))
	relative <- as.logical(rep_len(relative, ncol(table$table)))
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_kdtree_range_search, query, table, tolerance,
		relative, referent, as.integer(num.threads))
}

