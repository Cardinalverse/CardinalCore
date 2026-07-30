
#### Compute differences
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
