
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
	if ( !is.null(y) && (nrow(x) != nrow(y) ) )
		stop("nrow(x) [", nrow(x), "] must equal ",
			"nrow(y) [", nrow(y), "]")
	if ( !is.null(group) ) {
		if ( reorder ) {
			group <- as.integer(match(group, sort(unique(group))) - 1L)
		} else {
			group <- as.integer(match(group, unique(group)) - 1L)
		}
	}
	if ( !is.null(weights) )
		weights <- matrix(as.double(weights), nrow=nrow(x))
	.Call(C_do_col_sums, x, as.integer(num.threads))
}
