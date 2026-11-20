\name{kdsearch}

\alias{kdtree}
\alias{kdsearch}
\alias{kdsearch_agg}
\alias{kdsearch_first}
\alias{kdsearch_last}
\alias{knnsearch}

\title{K-Dimensional Search}

\description{
	Use a kd-tree for range searches and nearest neighbor searches on a table and aggregate over the matches.
}

\usage{
# K-D tree
kdtree(table)

# Range search
kdsearch(
	query,
	table,
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	nomatch = NA_integer_,
	split = TRUE,
	num.threads = 0L)

kdsearch_agg(
	query,
	table,
	values,
	stat = c("sum", "prod", "max", "min", "mean", "var"),
	tolerance = 0,
	relative = !missing(relative_to),
	relative_to = c("query", "table"),
	num.threads = 0L)

kdsearch_first(query, table, \dots)

kdsearch_last(query, table, \dots)

# Nearest neighbor search
knnsearch(
	query,
	table,
	k = 1L,
	metric = c("Euclidean", "Manhattan", "Maximum"),
	split = TRUE,
	num.threads = 0L)

}

\arguments{
	\item{query}{A numeric matrix or data frame with rows as queries.}

	\item{table}{A \code{kdtree} returned by \code{kdtree()}, or a numeric matrix or data frame with rows as records to be matched against.}

	\item{values}{A numeric vector of values to be aggregated, with length equal to length of \code{table}.}

	\item{stat}{The summary statistic for aggregation.}

	\item{tolerance}{A numeric vector tolerance for matching. Must all be >= 0.}

	\item{relative}{Should relative differences be used?}
	j
	\item{relative_to}{The reference for relative differences.}

	\item{nomatch}{The value to be returned when no matches are found, coerced to an integer.}

	\item{k}{The number of nearest neighbors to find for each row in \code{query}.}

	\item{metric}{Distance metric to use for nearest neighbors.}

	\item{split}{Should the results be split into a list with the same length as the number of rows in \code{query}?}

	\item{num.threads}{The number of threads, parallelized over \code{query}, where 0 means no thread dispatch at all, and 1 or greater dispatches work to that number of worker threads.}

	\item{\dots}{Additional arguments to pass to \code{kdsearch_agg}.}
}

\details{
    A kd-tree is a multidimensional generalization of a binary search tree. Building the search tree is O(n * log n) and searching for a single query is O(log n). If multiple searches or aggregations will be applied to the same data, it can be much faster to build the tree once with \code{kdtree()}.

    These functions are implemented in C++ with multithreading support.
}

\value{
	For \code{kdsearch()}, defaults to a list with the same length as rows in \code{query}, where entries as vectors of indices of matches.

	For \code{knnsearch()}, defaults to a list with following components:

	\item{index}{A list with the same length as rows in \code{query}, where entries are vectors of indices of nearest neighbors.}
	\item{distance}{A list with the same length as rows in \code{query}, where entries are vectors of distances to nearest neighbors.}
	
	For \code{kdtree()}, a \code{kdtree} object for more efficient searches.
}

\author{Kylie A. Bemis}

\seealso{
	\code{\link{bsearch}}
}

\examples{
table <- expand.grid(x=1:10, y=1:10)
q <- rbind(c(1.11, 2.22), c(3.33, 4.44))

knnsearch(q, table, k=3)
}

\keyword{tree}
\keyword{spatial}
\keyword{utilities}
