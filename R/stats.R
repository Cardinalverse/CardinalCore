
#### Streaming means
## ------------------

stream_means <- function(means, nobs)
{
	structure(means, nobs=nobs, class="stream_means")
}

stream_vars <- function(vars, means, nobs)
{
	structure(vars, means=means, nobs=nobs, class="stream_vars")
}

merge_means <- function(x, y, ...)
{
	if ( !inherits(x, "stream_means") || !inherits(y, "stream_means") )
		stop("x and y must be stream_means object")
	.Call(C_do_stream_means, x, y)
}

merge_vars <- function(x, y, ...)
{
	if ( !inherits(x, "stream_vars") || !inherits(y, "stream_vars") )
		stop("x and y must be stream_vars object")
	.Call(C_do_stream_vars, x, y)
}

#### Compute column sums
## ---------------------

col_sums <- function(
	x,
	group = NULL,
	reorder = TRUE,
	num.threads = 1)
{
	if ( !is.null(group) ) {
		if ( length(group) != NROW(x) )
			stop("length(group) [", length(group), "] and ",
				"NROW(x) [", NROW(x), "] must be equal")
		if ( anyNA(group) )
			stop("missing values for 'group'")
		if ( reorder ) {
			ugroup <- sort(unique(group))
		} else {
			ugroup <- unique(group)
		}
		group <- as.integer(match(group, ugroup) - 1L)
		ngroups <- length(ugroup)
	}
	if ( is.null(group) ) {
		.Call(C_do_col_sums, x, as.integer(num.threads))
	} else {
		ans <- .Call(C_do_col_sums, x, 
			group, ngroups, as.integer(num.threads))
		rownames(ans) <- as.character(ugroup)
		ans
	}
}
