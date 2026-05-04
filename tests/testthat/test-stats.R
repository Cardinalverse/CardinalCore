require(CardinalCore)
require(testthat)

test_that("col_sums works", {

	set.seed(1)
	nr <- 10000
	nc <- 5000
	x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)

	expect_equal(col_sums(x), colSums(x, na.rm=TRUE))

	matter::mem(x)
	bench::mark(col_sums(x, 1), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, 2), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, 4), colSums(x, na.rm=TRUE))

})
