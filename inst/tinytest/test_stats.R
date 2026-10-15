require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# stream_stats
x <- stream_stats(runif(10), stat="mean")
y <- stream_stats(runif(10), stat="mean")
z <- merge_stats(x, y)
group_stats(x, rep(1:2, 5))
group_stats(x, rep(1:5, 2))
merge_stats(x, z)

# Setup example matrix
set.seed(1)
nr <- 1e6
nc <- 128
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
g <- sample(8, nrow(x), replace=TRUE)

c2 <- function(x) {
	attributes(x) <- attributes(x)[c("dim", "dimnames")]
	x
}

expect_equal(c(col_stats(x, "sum")), colSums(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "sum")), colSums(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "sum", num.threads=2)), colSums(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "sum", num.threads=4)), colSums(x, na.rm=TRUE))

expect_equal(c(col_stats(x, "mean")), colMeans(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "mean")), colMeans(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "mean", num.threads=2)), colMeans(x, na.rm=TRUE))
bench::mark(c(col_stats(x, "mean", num.threads=4)), colMeans(x, na.rm=TRUE))

expect_equal(c(col_stats(x, "var")), apply(x, 2, var, na.rm=TRUE))
bench::mark(c(col_stats(x, "var")), apply(x, 2, var, na.rm=TRUE))
bench::mark(c(col_stats(x, "var", num.threads=2)), apply(x, 2, var, na.rm=TRUE))
bench::mark(c(col_stats(x, "var", num.threads=4)), apply(x, 2, var, na.rm=TRUE))

expect_equal(c2(col_stats(x, "sum", group=g)), rowsum(x, group=g))
bench::mark(c2(col_stats(x, "sum", group=g, num.threads=1)), rowsum(x, group=g))
bench::mark(c2(col_stats(x, "sum", group=g, num.threads=2)), rowsum(x, group=g))
bench::mark(c2(col_stats(x, "sum", group=g, num.threads=4)), rowsum(x, group=g))

xt <- t(x)
expect_equal(c(row_stats(xt, "sum")), rowSums(xt, na.rm=TRUE))
bench::mark(c(row_stats(xt, "sum")), rowSums(xt, na.rm=TRUE))
bench::mark(c(row_stats(xt, "sum", num.threads=2)), rowSums(xt, na.rm=TRUE))
bench::mark(c(row_stats(xt, "sum", num.threads=4)), rowSums(xt, na.rm=TRUE))

xt <- t(x)
expect_equal(c2(row_stats(xt, "sum", group=g)), rowsum(x, group=g))
bench::mark(c2(row_stats(xt, "sum", group=g, num.threads=1)), rowsum(x, group=g))
bench::mark(c2(row_stats(xt, "sum", group=g, num.threads=2)), rowsum(x, group=g))
bench::mark(c2(row_stats(xt, "sum", group=g, num.threads=4)), rowsum(x, group=g))

# test expression templates
test <- function(x, i) {
	.Call(CardinalCore:::C_do_test_expression, x, i)
}
x <- runif(1e6)
i <- sample(1e3)
i0 <- i - 1L
expect_equal(log1p(x + x)[i], test(x, i0))

bench::mark(log1p(x + x)[i], test(x, i0))

