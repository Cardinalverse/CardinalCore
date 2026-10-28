require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# col_stats and row_stats
set.seed(1)
nr <- 1e6
nc <- 32L
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
x[sample(100)] <- NA_real_
gr <- sample(6L, nr, replace=TRUE)
gc <- sample(6L, nc, replace=TRUE)
bench::mark(
	col_stats(x, "sum", num.threads=1),
	col_stats(x, "sum", num.threads=2),
	col_stats(x, "sum", num.threads=4))
bench::mark(
	col_stats(x, "sum", gr, num.threads=1),
	col_stats(x, "sum", gr, num.threads=2),
	col_stats(x, "sum", gr, num.threads=4))

