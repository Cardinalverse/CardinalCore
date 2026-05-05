require(CardinalCore)
require(testthat)

test_that("col_sums works", {

	set.seed(1)
	nr <- 10000
	nc <- 5000
	x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
	g <- sample(8, nrow(x), replace=TRUE)

	expect_equal(col_sums(x), colSums(x, na.rm=TRUE))
	expect_equal(col_sums(x, group=g), rowsum(x, group=g))

	bench::mark(col_sums(x, num.threads=1), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=2), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=4), colSums(x, na.rm=TRUE))

	bench::mark(
		col_sums(x, group=g, num.threads=1), 
		rowsum(x, group=g, na.rm=TRUE))
	bench::mark(
		col_sums(x, group=g, num.threads=2), 
		rowsum(x, group=g, na.rm=TRUE))
	bench::mark(
		col_sums(x, group=g, num.threads=4), 
		rowsum(x, group=g, na.rm=TRUE))

})
