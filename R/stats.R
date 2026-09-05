
#### Streaming means
## ------------------

stream_stats <- function(data = NA_real_, 
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	nobs = rep_len(1L, length(data)), ...)
{
	stat <- match.arg(stat)
	if ( stat == "var" ) {
		mean <- data
		data <- rep(NA_real_, length(data))
		structure(data, mean=mean, nobs=nobs, stat=stat, class="stream_stats")
	} else {
		structure(data, nobs=nobs, stat=stat, class="stream_stats")
	}
}

merge_stats <- function(x, y)
{
	if ( !inherits(x, "stream_stats") || !inherits(y, "stream_stats") )
		stop("'x' and 'y' must be stream_stats object")
	if ( length(x) != length(y) )
		stop("length(x) [", length(x), "] and ",
			"length(y) [", length(y), "] must be equal")
	if ( attr(x, "stat") != attr(y, "stat") )
		stop("attr(x, 'stat') must match attr(y, 'stat')")
	if ( attr(x, "stat") == "mean" ) {
		.Call(C_do_merge_means, x, y)
	} else if ( attr(x, "stat") == "var" ) {
		.Call(C_do_merge_vars, x, y)
	} else {
		.Call(C_do_merge_stats, x, y)
	}
}

group_stats <- function(x, group, reorder = TRUE)
{
	if ( !inherits(x, "stream_stats") )
		stop("'x' must be stream_stats object")
	if ( length(group) != length(x) )
		stop("length(group) [", length(group), "] and ",
			"length(x) [", length(x), "] must be equal")
	if ( anyNA(group) )
		stop("missing values in 'group'")
	if ( reorder ) {
		ugroup <- sort(unique(group))
	} else {
		ugroup <- unique(group)
	}
	group <- as.integer(match(group, ugroup) - 1L)
	if ( attr(x, "stat") == "mean" ) {
		.Call(C_do_group_means, x, group, length(ugroup))
	} else if ( attr(x, "stat") == "var" ) {
		.Call(C_do_group_vars, x, group, length(ugroup))
	} else {
		.Call(C_do_group_stats, x, group, length(ugroup))
	}
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
			stop("missing values in 'group'")
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
