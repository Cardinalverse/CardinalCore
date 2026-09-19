
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
	if ( is.unsorted(table) )
		stop("'table' must be sorted")
	if ( is.double(query) && is.integer(table) )
		table <- as.double(table)
	if ( is.integer(query) && is.double(table) )
		query <- as.double(query)
	relative <- isTRUE(relative)
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_bsearch, query, table,
		as.double(tolerance), relative, referent, as.integer(nomatch))
}

bsearch_agg <- function(
	query,
	table,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"))
{
	stat <- match.arg(stat)
	if ( is.unsorted(table) )
		stop("'table' must be sorted")
	if ( length(table) != length(values) )
		stop("length of 'values' must equal length of 'table'")
	if ( is.double(query) && is.integer(table) )
		table <- as.double(table)
	if ( is.integer(query) && is.double(table) )
		query <- as.double(query)
	values <- as.double(values)
	relative <- isTRUE(relative)
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_bsearch_aggregate, query, table, values, stat,
		as.double(tolerance), relative, referent)
}

bsearch_first <- function(query, table, ...)
{
	bsearch_agg(query, table, values=seq_along(table), stat="min", ...)
}

bsearch_last <- function(query, table, ...)
{
	bsearch_agg(query, table, values=seq_along(table), stat="max", ...)
}

#### Binary multisearch
## ---------------------

msearch <- function(
	query,
	tables,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_)
{
	if ( is.null(dim(query)) ) {
		query <- t(query)
	} else {
		query <- as.matrix(query)
	}
	if ( ncol(query) != length(tables) )
		stop("number of columns in 'query' must match length of 'tables'")
	tolerance <- rep_len(tolerance, ncol(query))
	relative <- rep_len(relative, ncol(query))
	hits <- matrix(nomatch, nrow=nrow(query), ncol=ncol(query))
	dimnames(hits) <- dimnames(query)
	for ( i in seq_along(tables) )
		hits[,i] <- bsearch(
			query[,i],
			tables[[i]],
			tolerance=tolerance[i],
			relative=relative[i],
			relative_to=relative_to,
			nomatch=nomatch)
	hits
}

msearch_agg <- function(
	query,
	tables,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"))
{
	if ( is.null(dim(query)) ) {
		query <- t(query)
	} else {
		query <- as.matrix(query)
	}
	if ( ncol(query) != length(tables) )
		stop("number of columns in query must match length of tables")
	if ( !identical(lengths(tables), lengths(values)) )
		stop("lengths of 'values' must equal lengths of 'tables'")
	tolerance <- rep_len(tolerance, ncol(query))
	relative <- rep_len(relative, ncol(query))
	aggs <- matrix(nomatch, nrow=nrow(query), ncol=ncol(query))
	dimnames(aggs) <- dimnames(query)
	for ( i in seq_along(tables) )
		aggs[,i] <- bsearch_agg(
			query,
			tables[[i]],
			values[[i]],
			stat=stat,
			tolerance=tolerance[i],
			relative=relative[i],
			relative_to=relative_to)
	aggs
}

msearch_first <- function(query, tables, ...)
{
	msearch_agg(query, tables, values=lapply(tables, seq_along), stat="min", ...)
}

msearch_last <- function(query, tables, ...)
{
	msearch_agg(query, tables, values=lapply(tables, seq_along), stat="max", ...)
}

#### Kd-tree search
## ----------------

kdtree <- function(table)
{
	if ( inherits(table, "kdtree") )
		return(table)
	.Call(C_do_kdtree_build, as.matrix(table))
}

kdsearch <- function(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
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
	if ( ncol(query) != ncol(table$table) )
		stop("'query' must have the same number of columns as 'table'")
	if ( anyNA(tolerance) )
		stop("'tolerance' must not contain NAs")
	if ( anyNA(relative) )
		stop("'relative' must not contain NAs")
	if ( is.integer(query) && is.double(table$table) )
		storage.mode(query) <- "double"
	if ( is.double(query) && is.integer(table$table) )
		storage.mode(table$table) <- "double"
	tolerance <- as.double(rep_len(tolerance, ncol(table$table)))
	relative <- as.logical(rep_len(relative, ncol(table$table)))
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_kdtree_range_search, query, table,
		tolerance, relative, referent, as.integer(num.threads))
}

kdsearch_agg <- function(
	query,
	table,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	num.threads = 1)
{
	stat <- match.arg(stat)
	if ( !inherits(table, "kdtree") )
		table <- kdtree(table)
	if ( is.null(dim(query)) ) {
		query <- t(query)
	} else {
		query <- as.matrix(query)
	}
	if ( length(values) != nrow(table$table) )
		stop("length of 'values' must equal number of rows in 'table'")
	if ( ncol(query) != ncol(table$table) )
		stop("'query' must have the same number of columns as 'table'")
	if ( anyNA(tolerance) )
		stop("'tolerance' must not contain NAs")
	if ( anyNA(relative) )
		stop("'relative' must not contain NAs")
	if ( is.integer(query) && is.double(table$table) )
		storage.mode(query) <- "double"
	if ( is.double(query) && is.integer(table$table) )
		storage.mode(table$table) <- "double"
	values <- as.double(values)
	tolerance <- as.double(rep_len(tolerance, ncol(table$table)))
	relative <- as.logical(rep_len(relative, ncol(table$table)))
	referent <- c("query"=0L, "table"=1L)[match.arg(relative_to)]
	.Call(C_do_kdtree_range_aggregate, query, table, values, stat,
		tolerance, relative, referent, as.integer(num.threads))
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

