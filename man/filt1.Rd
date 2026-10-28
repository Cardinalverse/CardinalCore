\name{filt1}

\alias{filt1}

\alias{filt1_mean}
\alias{filt1_conv}

\title{Smoothing Filters in 1D}

\description{
    Smooth a uniformly sampled 1D signal.
}

\usage{
# Moving average filter
filt1_mean(y, k = 5L)

# Linear convolution filter
filt1_conv(y, w = rep(1, 5L))
}

\arguments{
	\item{y}{A numeric vector.}

    \item{k}{The number of samples in the smoothing window. Must be positive. Must be odd.}

    \item{w}{The weights of the linear convolution kernel. Length must be odd.}
}

\details{
    All filters documented here extend the signal by replicating the endpoints.

    \code{filt1_mean()} performs mean filtering in O(n) time. This is fast and especially useful for calculating other filters that can be constructed as a combination of mean filters.

    \code{filt1_conv()} performs a convolution filtering with arbitary weights. The weights are internally normalized to sum to 1.
}

\value{
    A numeric vector the same length as \code{y} with the smoothed result.
}

\author{Kylie A. Bemis}

\examples{
set.seed(1)
x <- seq(from=0, to=6 * pi, length.out=5000)
y <- sin(x) + 0.6 * sin(2.6 * t)
y <- y + runif(length(y))
z <- filt1_mean(y, k=25)

plot(y, type="l")
lines(z, col="red")
}

\keyword{ts}
\keyword{smooth}
