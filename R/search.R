
#### Binary search
## ----------------

bsearch <- function(query, x, tolerance = 0,
	relative = FALSE, nearest = FALSE, nomatch = NA_integer_)
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

