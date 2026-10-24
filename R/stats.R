
#### Streaming statistics
## ----------------------

stream_stats <- function(x = numeric(), nobs = 1L, 
	stat = c("sum", "prod", "max", "min", "mean", "var"))
{
	stat <- match.arg(stat)
	nobs = rep_len(nobs, length(x))
	if ( stat == "var" ) {
		if ( any(nobs > 1L) )
			stop("initializing nobs > 1 not supported for stat 'var'")
		structure(rep(NA_real_, length(x)), mean=as.double(x), 
			nobs=nobs, stat=stat, class="stream_stats")
	} else {
		structure(as.double(x), 
			nobs=nobs, stat=stat, class="stream_stats")
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
		cat("nobs():", n, "\n")
	}
}

nobs.stream_stats <- function(object, ...) {
	attr(object, "nobs")
}

merge_stats <- function(x, y, ...)
{
	if ( !inherits(x, "stream_stats") )
		x <- stream_stats(x, ...)
	if ( !inherits(y, "stream_stats") )
		y <- stream_stats(y, ...)
	if ( length(x) != length(y) )
		stop("length(x) [", length(x), "] and ",
			"length(y) [", length(y), "] must be equal")
	if ( attr(x, "stat") != attr(y, "stat") )
		stop("'x' is a ", attr(x, "stat"), " but 'y' is a ", attr(y, "stat"))
	.Call(C_do_merge_stats, x, y)
}

group_stats <- function(x, group, reorder = TRUE, ...)
{
	if ( !inherits(x, "stream_stats") )
		x <- stream_stats(x, ...)
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
	.Call(C_do_group_stats, x, group, ugroup)
}

collapse_by <- function(x, by, stat, count = "nobs")
{
	if ( by %notin% colnames(x) )
		stop("missing 'by' column [", by, "]")
	group <- x[,by,drop=TRUE]
	if ( count %in% colnames(x) ) {
		nobs <- x[,count,drop=TRUE]
	} else {
		nobs <- as.integer(!is.na(group))
	}
	options <- c("sum", "prod", "max", "min", "mean", "var")
	stat <- resolve_named_opts(stat, options, names=colnames(x), 
		exclude=c(by, count), duplicates.ok=FALSE)
	out <- setNames(vector("list", length(stat)), names(stat))
	for ( i in seq_along(out) ) {
		nm <- names(stat)[i]
		out[[nm]] <- group_stats(x[,nm,drop=TRUE], 
			group=group, stat=stat[i], nobs=nobs)
	}
	out[[count]] <- group_stats(nobs, group=group, stat="sum")
	out[[by]] <- attr(out[[count]], "group")
	out
}

resolve_named_opts <- function(arg, options, names, 
	exclude = character(), duplicates.ok = FALSE)
{
	i <- pmatch(arg, options, duplicates.ok=TRUE)
	unsupported <- which(is.na(i))
	if ( length(unsupported) > 0 ) {
		badoptions <- unique(arg[unsupported])
		stop("unsupported option(s) [", paste0(badoptions, collapse=", "), "]")
	}
	if ( is.null(names(arg)) )
		names(arg) <- rep.int("", length(arg))
	default <- which(names(arg) == "")
	if ( length(default) > 1 ) {
		stop("only one option can be unnamed (default)")
	} else if ( length(default) > 0 ) {
		cols <- setdiff(names, c(names(arg), exclude))
		cols <- setNames(rep.int(arg[default], length(cols)), cols)
		arg <- c(arg[-default], cols)
	}
	names(arg) <- match.arg(names(arg), names, several.ok=TRUE)
	if ( anyDuplicated(names(arg)) ) {
		dups <- names(arg)[duplicated(names(arg))]
		stop("duplicate named option(s) [", paste0(dups, collapse=", "), "]")
	}
	arg
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

