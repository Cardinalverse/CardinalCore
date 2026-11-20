\name{qorder}

\alias{qorder}
\alias{qselect}
\alias{qmedian}
\alias{qmad}

\title{Internal Sorting}

\description{
    Internal utilities for sorting and selection.
}

\usage{
# QuickSelect-based k-th order statistic
qselect(x, k)

# QuickSort-based index ordering
qorder(x)

# Median
qmedian(x)

# Median-absolute deviation
qmad(x, center = qmedian(x), constant = 1.4826)
}

\arguments{
	\item{x}{A numeric vector.}

    \item{k}{The k-th order statistics.}

    \item{center}{The center for computing absolute deviations.}

    \item{constant}{The scaling factor.}
}

\details{
    These are internal utilities intended for testing C++ sorting and selection used in native implementations of other functions.

    Sort order is always ascending and \code{NA} values always sort last.
}

\value{
    For \code{qselect} and \code{qorder}, an integer vector giving the sort indices or the indices of the k-th order statistics.

    For \code{qmedian} and \code{qmad}, a numeric scalar.
}

\author{Kylie A. Bemis}

\examples{
set.seed(1)
x <- runif(10)

all.equal(qorder(x), order(x))
all.equal(qmedian(x), median(x))
all.equal(qmad(x), mad(x))
}

\keyword{utilities}

