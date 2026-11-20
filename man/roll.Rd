\name{roll}

\alias{roll}
\alias{roll_apply}

\title{Rolling Apply}

\description{
    A simple rolling apply over vectors.
}

\usage{
# Rolling windows
roll(y, k = 5L)

# Rolling function application
roll_apply(y, k = 5L, FUN = NULL, \dots, simplify = TRUE)
}

\arguments{
	\item{y}{A numeric vector.}

    \item{k}{The number of samples in the rolling window. Must be positive. Must be odd.}

    \item{FUN}{The function to apply to the windows.}

    \item{\dots}{Arguments to pass to \code{FUN}.}

    \item{simplify}{Should \code{\link{simplify2array}} be applied to the output?}
}

\details{
    Endpoints are repeated for padding. These implementations make no attempt at efficiency. They are provided primarily for convenience and testing purposes.
}

\value{
    For \code{roll}, a list of vectors (each of length \code{k}).

    For \code{roll_apply}, a vector with elements corresponding to the result of \code{FUN} applied to the rolling windows.
}

\author{Kylie A. Bemis}

\examples{
set.seed(1)
x <- runif(10)

roll(x, 3)
roll_apply(x, 3, mean)
}

\keyword{utilities}
