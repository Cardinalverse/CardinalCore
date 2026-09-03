require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# Setup example streams
nr <- 100
nc <- 10
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
y <- matrix(runif(nr * nc), nrow=nr, ncol=nc)

sx <- apply(x, 2L, var)
mx <- colMeans(x)
nx <- rep(nr, nc)
sy <- apply(y, 2L, var)
my <- colMeans(y)
ny <- rep(nr, nc)

smx <- stream_means(mx, nx)
smy <- stream_means(my, ny)
ssx <- stream_vars(sx, mx, nx)
ssy <- stream_vars(sy, my, ny)

merge_means(smx, smy)
colMeans(rbind(x, y))

merge_vars(ssx, ssy)
apply(rbind(x, y), 2L, var)

# Setup example matrix
set.seed(1)
nr <- 1e6
nc <- 72
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
g <- sample(8, nrow(x), replace=TRUE)

# col_sums
expect_equal(col_sums(x), colSums(x, na.rm=TRUE))

bench::mark(col_sums(x, num.threads=1), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=3), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=6), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=12), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=24), colSums(x, na.rm=TRUE))

# test expression templates
test <- function(x, i) {
	.Call(CardinalCore:::C_do_test_expression, x, i)
}
x <- runif(1e6)
i <- sample(1e3)
i0 <- i - 1L
expect_equal(log1p(x + x)[i], test(x, i0))

bench::mark(log1p(x + x)[i], test(x, i0))

