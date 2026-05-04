
#### Compute differences
## ---------------------

col_sums <- function(
	x,
	group,
	center = FALSE,
	scale = FALSE,
	num.threads = 1)
{
	.Call(C_do_col_sums, x, as.integer(num.threads))
}
