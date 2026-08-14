require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# Setup example matrix
set.seed(1)
nr <- 1e6
nc <- 64
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
g <- sample(8, nrow(x), replace=TRUE)

# col_sums
expect_equal(col_sums(x), colSums(x, na.rm=TRUE))

bench::mark(col_sums(x, num.threads=1), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=2), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=4), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=8), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=16), colSums(x, na.rm=TRUE))
bench::mark(col_sums(x, num.threads=32), colSums(x, na.rm=TRUE))

# test expression templates
test <- function(x, i) {
	.Call(CardinalCore:::C_do_test_expression, x, i)
}
x <- runif(1e6)
i <- sample(1e3)
i0 <- i - 1L
expect_equal(log1p(x + x)[i], test(x, i0))

bench::mark(log1p(x + x)[i], test(x, i0))

