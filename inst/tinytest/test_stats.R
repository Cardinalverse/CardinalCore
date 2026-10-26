require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# stream_stats methods
x <- stream_stats(runif(10L), stat="mean")
y <- stream_stats(runif(10L), stat="var")
dim(x) <- c(5,2)
dim(y) <- c(5,2)
expect_equal(dim(x), c(5,2))
expect_equal(dim(y), c(5,2))
expect_equal(dim(nobs(x)), c(5,2))
expect_equal(dim(nobs(y)), c(5,2))
dim(x) <- NULL
dim(y) <- NULL
expect_null(dim(x))
expect_null(dim(y))
expect_inherits(t(x), "stream_stats")
expect_inherits(x[1:5], "stream_stats")
x[1:5] <- 1
y[11:15] <- 1
expect_inherits(x, "stream_stats")
expect_inherits(y, "stream_stats")
expect_length(x, 10L)
expect_length(y, 15L)

# merge_stats
set.seed(1)
x <- runif(10)
y <- runif(10)
expect_equal(c(merge_stats(x, y, stat="sum")), rowSums(cbind(x, y)))
expect_equal(c(merge_stats(x, y, stat="mean")), rowMeans(cbind(x, y)))
expect_equal(c(merge_stats(x, y, stat="var")), apply(cbind(x, y), 1L, var))

# group_stats
set.seed(1)
n <- 100L
x <- runif(n)
y <- runif(n)
gx <- sample(6L, n, replace=TRUE)
gy <- sample(6L, n, replace=TRUE)
s1x <- group_stats(x, group=gx, stat="mean")
s1y <- group_stats(y, group=gy, stat="mean")
expect_equal(c(s1x), c(tapply(x, gx, mean)))
expect_equal(c(s1y), c(tapply(y, gy, mean)))
expect_equal(c(merge_stats(s1x, s1y)), c(tapply(c(x, y), c(gx, gy), mean)))
s2x <- group_stats(x, group=gx, stat="var")
s2y <- group_stats(y, group=gy, stat="var")
expect_equal(c(s2x), c(tapply(x, gx, var)))
expect_equal(c(s2y), c(tapply(y, gy, var)))
expect_equal(c(merge_stats(s2x, s2y)), c(tapply(c(x, y), c(gx, gy), var)))

mtcars
as.data.frame(collapse_by(mtcars, "gear", stat="mean"))
as.data.frame(collapse_by(mtcars, "gear", stat=c(disp="mean", "max")))
collapse_by(mtcars, "gear", stat=c(disp="mean", mpg="mean"))

l <- collapse_by(mtcars, "gear", stat=c(disp="mean", mpg="mean"))
d <- as.data.frame(collapse_by(mtcars, "gear", stat=c(disp="mean", mpg="mean")))
# collapse_by(mtcars, "gear", stat=c(disp="mean", disp="mean"))

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

