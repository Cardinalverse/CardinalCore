
#### Streaming statistics
## ----------------------

stream_stats <- function(x = numeric(), 
	stat = c("sum", "prod", "max", "min", "mean", "var"))
{
	stat <- match.arg(stat)
	nobs = rep_len(1L, length(x))
	if ( stat == "var" ) {
		mean <- x
		x <- rep(NA_real_, length(x))
		structure(x, mean=mean, nobs=nobs, stat=stat, class="stream_stats")
	} else {
		structure(x, nobs=nobs, stat=stat, class="stream_stats")
	}
}

print.stream_stats <- function(x, nobs.max = 6L, ...)
{
	cat("Streaming ", attr(x, "stat"), "(s)\n", sep="")
	n <- attr(x, "nobs")
	s <- c(x)
	attributes(s) <- list(dim=dim(x), dimnames=dimnames(x), names=names(x))
	print.default(s, ...)
	if ( length(n) > nobs.max ) {
		cat("nobs():", head(as.vector(n), nobs.max), " ...\n")
	} else {
		cat("nobs():", n)
	}
}

nobs.stream_stats <- function(object, ...) {
	attr(object, "nobs")
}

merge_stats <- function(x, y)
{
	if ( !inherits(x, "stream_stats") || !inherits(y, "stream_stats") )
		stop("'x' and 'y' must be stream_stats object")
	if ( length(x) != length(y) )
		stop("length(x) [", length(x), "] and ",
			"length(y) [", length(y), "] must be equal")
	if ( attr(x, "stat") != attr(y, "stat") )
		stop("'x' is a ", attr(x, "stat"), " but 'y' is a ", attr(y, "stat"))
	.Call(C_do_merge_stats, x, y)
}

group_stats <- function(x, group, reorder = TRUE)
{
	if ( !inherits(x, "stream_stats") )
		stop("'x' must be a stream_stats object")
	if ( length(group) != length(x) )
		stop("length(group) [", length(group), "] and ",
			"length(x) [", length(x), "] must be equal")
	if ( anyNA(group) )
		stop("missing values in 'group'")
	if ( is.factor(group) ) {
		ugroup <- levels(group)
		group <- as.integer(group) - 1L
	} else {
		ugroup <- if (reorder) sort(unique(group)) else unique(group)
		group <- as.integer(match(group, ugroup) - 1L)
	}
	.Call(C_do_group_stats, x, group, as.character(ugroup))
}

#### Matrix statistics
## -------------------

col_stats <- function(
	x, 
	stat, 
	group = NULL,
	reorder = TRUE,
	num.threads = 1)
{
	stat <- match.arg(stat, c("sum", "prod", "max", "min", "mean", "var"))
	if ( !is.null(group) )
	{
		if ( length(group) != NROW(x) )
			stop("length(group) [", length(group), "] and ",
				"NROW(x) [", NROW(x), "] must be equal")
		if ( anyNA(group) )
			stop("missing values in 'group'")
		if ( is.factor(group) ) {
			ugroup <- levels(group)
			group <- as.integer(group) - 1L
		} else {
			ugroup <- if (reorder) sort(unique(group)) else unique(group)
			group <- as.integer(match(group, ugroup)) - 1L
		}
	}
	if ( is.null(group) ) {
		.Call(C_do_col_stats, as.matrix(x), stat, as.integer(num.threads))
	} else {
		.Call(C_do_col_group_stats, 
			as.matrix(x), stat, group, ugroup, as.integer(num.threads))
	}
}

row_stats <- function(
	x, 
	stat, 
	group = NULL,
	reorder = TRUE,
	num.threads = 1)
{
	stat <- match.arg(stat, c("sum", "prod", "max", "min", "mean", "var"))
	if ( !is.null(group) )
	{
		if ( length(group) != NCOL(x) )
			stop("length(group) [", length(group), "] and ",
				"NCOL(x) [", NCOL(x), "] must be equal")
		if ( anyNA(group) )
			stop("missing values in 'group'")
		if ( is.factor(group) ) {
			ugroup <- levels(group)
			group <- as.integer(group) - 1L
		} else {
			ugroup <- if (reorder) sort(unique(group)) else unique(group)
			group <- as.integer(match(group, ugroup)) - 1L
		}
	}
	if ( is.null(group) ) {
		.Call(C_do_row_stats, as.matrix(x), stat, as.integer(num.threads))
	} else {
		.Call(C_do_row_group_stats, 
			as.matrix(x), stat, group, ugroup, as.integer(num.threads))
	}
}

