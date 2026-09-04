
#### Binary search
## ----------------

bsearch <- function(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	which = c("nearest", "all"),
	nomatch = NA_integer_)
{
	if ( is.unsorted(table) )
		stop("'table' must be sorted")
	if ( is.double(query) && is.integer(table) )
		table <- as.double(table)
	if ( is.integer(query) && is.double(table) )
		query <- as.double(query)
	relative <- isTRUE(relative)
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	which <- match.arg(which)
	if ( which == "nearest" ) {
		.Call(C_do_bsearch, query, table, as.double(tolerance),
			relative, referent, as.integer(nomatch))
	} else {
		.Call(C_do_rsearch, query, table, as.double(tolerance),
			relative, referent, as.integer(nomatch))
	}
}

kdtree <- function(table)
{
	if ( inherits(table, "kdtree") )
		return(table)
	if ( is.null(dim(table)) )
		table <- as.matrix(table)
	.Call(C_do_kdtree_build, table)
}

kdsearch <- function(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	which = c("all", "first", "last"),
	nomatch = NA_integer_,
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
	which <- match.arg(which)
	if ( which == "all" ) {
		.Call(C_do_kdtree_range_search, query, table, tolerance,
			relative, referent, as.integer(num.threads))
	} else if ( which == "first" ) {
		.Call(C_do_kdtree_range_find_first, query, table, tolerance,
			relative, referent, as.integer(nomatch), as.integer(num.threads))
	} else if ( which == "last" ) {
		.Call(C_do_kdtree_range_find_last, query, table, tolerance,
			relative, referent, as.integer(nomatch), as.integer(num.threads))
	}
}

knnsearch <- function(
	query,
	table,
	k = 1L,
	metric = c("Euclidean", "Manhattan", "Maximum"),
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
	k <- rep_len(as.integer(k), nrow(query))
	metric <- match.arg(metric)
	metric <- c("Manhattan"=0L, "Euclidean"=1L, "Maximum"=2L)[metric]
	.Call(C_do_kdtree_knn_search, query, table, k,
		metric, as.integer(num.threads))
}

