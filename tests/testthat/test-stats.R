require(CardinalCore)
require(testthat)

test_that("col_sums works", {

	set.seed(1)
	nr <- 1e6
	nc <- 32
	x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
	g <- sample(8, nrow(x), replace=TRUE)

	expect_equal(col_sums(x), colSums(x, na.rm=TRUE))
	# expect_equal(col_sums(x, group=g), rowsum(x, group=g))

	bench::mark(col_sums(x, num.threads=1), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=2), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=4), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=8), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=16), colSums(x, na.rm=TRUE))
	bench::mark(col_sums(x, num.threads=32), colSums(x, na.rm=TRUE))

	# bench::mark(
	# 	col_sums(x, group=g, num.threads=1), 
	# 	rowsum(x, group=g, na.rm=TRUE))
	# bench::mark(
	# 	col_sums(x, group=g, num.threads=2), 
	# 	rowsum(x, group=g, na.rm=TRUE))
	# bench::mark(
	# 	col_sums(x, group=g, num.threads=4), 
	# 	rowsum(x, group=g, na.rm=TRUE))
	#
})

test_that("test_expression", {

	test <- function(x, i) {
		.Call(CardinalCore:::C_do_test_expression, x, i)
	}
	x <- runif(100)
	i <- c(1L, 2L, 3L, 4L, 5L)

	log1p(x + x)[i]
	test(x, i - 1L)

	bench::mark(
		log1p(x + x)[i],
		test(x, i - 1L))

})


