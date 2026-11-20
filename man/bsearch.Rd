\name{bsearch}

\alias{bsearch}
\alias{bsearch_agg}
\alias{bsearch_first}
\alias{bsearch_last}

\alias{msearch}
\alias{msearch_agg}
\alias{msearch_first}
\alias{msearch_last}

\title{Binary Search}

\description{
    Use binary search to perform approximate matching and aggregation within a specified tolerance.
}

\usage{
# Binary search
bsearch(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_,
	num.threads = 0L)

bsearch(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_,
	num.threads = 0L)

bsearch_agg(
	query,
	table,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	num.threads = 0L)

bsearch_first(query, table, \dots)

bsearch_last(query, table, \dots)

# Multicolumn binary search
msearch(
	query,
	tables,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_,
	num.threads = 0L)

msearch_agg(
	query,
	tables,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	num.threads = 0L)

msearch_first(query, tables, \dots)

msearch_last(query, tables, \dots)
}

\arguments{
	\item{query}{A numeric vector of query values, or a matrix or data frame with columns of such vectors for multicolumn queries.}

	\item{table, tables}{A sorted vector of values to be matched against, or a list of such vectors for multicolumn queries.}

	\item{values}{A numeric vector of values to be aggregated, with length equal to length of \code{table}, or a list of such vectors for multicolumn queries.}

	\item{stat}{The summary statistic for aggregation.}

	\item{tolerance}{A scalar tolerance for matching. Must be >= 0.}

	\item{relative}{Should relative differences be used?}

	\item{relative_to}{The reference for relative differences.}

	\item{nomatch}{The value to be returned when no matches are found, coerced to an integer.}

	\item{num.threads}{The number of threads, parallelized over \code{query}, where 0 means no thread dispatch at all, and 1 or greater dispatches work to that number of worker threads.}

	\item{\dots}{Additional arguments to pass to \code{bsearch_agg}.}
}

\details{
    These functions are implemented in C++ with multithreading support.
}

\value{
    A numeric vector of the same length as \code{query}, or a matrix with a row for every row in \code{query} for multicolumn queries.
}

\author{Kylie A. Bemis}

\seealso{
    \code{\link{kdsearch}}
}

\examples{
ref <- c(1.11, 2.22, 3.33, 5.0, 5.1)

bsearch(2.22, ref) # 2
bsearch(3.0, ref) # NA
bsearch(3.0, ref, tolerance=Inf) # 3
bsearch(3.0, ref, tolerance=0.1, relative_to="table") # 3
bsearch(3.0, ref, tolerance=0.1, relative_to="query") # NA
}

\keyword{arith}
\keyword{utilities}
