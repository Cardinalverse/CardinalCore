
#### Compute differences
## ---------------------

col_sums <- function(
	x, y = NULL,
	group = NULL,
	reorder = TRUE,
	weights = NULL,
	abs = FALSE,
	p = 1,
	num.threads = 1)
{
	if ( !is.null(y) && (NROW(x) != NROW(y) ) )
		stop("NROW(x) [", NROW(x), "] and ",
			"NROW(y) [", NROW(y), "] must be equal")
	if ( !is.null(group) ) {
		if ( length(group) != NROW(x) )
			stop("length(group) [", length(group), "] and ",
				"NROW(y) [", NROW(y), "] must be equal")
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
	if ( !is.null(weights) )
		weights <- matrix(as.double(weights), nrow=nrow(x))
	if ( is.null(group) ) {
		.Call(C_do_col_sums, x, as.integer(num.threads))
	} else {
		ans <- .Call(C_do_col_sums_grouped, x, group, ngroups, as.integer(num.threads))
		rownames(ans) <- as.character(ugroup)
		ans
	}
}
