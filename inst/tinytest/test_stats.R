require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# stream_stats
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
s1x <- group_stats(x, gx, stat="sum")
s1y <- group_stats(y, gy, stat="sum")
expect_equal(c(s1x), c(tapply(x, gx, sum)))
expect_equal(c(s1y), c(tapply(y, gy, sum)))
expect_equal(c(merge_stats(s1x, s1y)), c(tapply(c(x, y), c(gx, gy), sum)))
s2x <- group_stats(x, gx, stat="mean")
s2y <- group_stats(y, gy, stat="mean")
expect_equal(c(s2x), c(tapply(x, gx, mean)))
expect_equal(c(s2y), c(tapply(y, gy, mean)))
expect_equal(c(merge_stats(s2x, s2y)), c(tapply(c(x, y), c(gx, gy), mean)))
s3x <- group_stats(x, gx, stat="var")
s3y <- group_stats(y, gy, stat="var")
expect_equal(c(s3x), c(tapply(x, gx, var)))
expect_equal(c(s3y), c(tapply(y, gy, var)))
expect_equal(c(merge_stats(s3x, s3y)), c(tapply(c(x, y), c(gx, gy), var)))

# col_stats and row_stats
set.seed(1)
nr <- 1000L
nc <- 32L
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
gr <- sample(6L, nr, replace=TRUE)
gc <- sample(6L, nc, replace=TRUE)
f <- function(INDEX, FUN) function(x) tapply(x, INDEX, FUN)
expect_equal(c(col_stats(x, "sum")), colSums(x))
expect_equal(c(row_stats(x, "sum")), rowSums(x))
expect_equal(c(col_stats(x, "mean")), colMeans(x))
expect_equal(c(row_stats(x, "mean")), rowMeans(x))
expect_equal(c(col_stats(x, "var")), apply(x, 2L, var))
expect_equal(c(row_stats(x, "var")), apply(x, 1L, var))
expect_equal(
	cbind(col_stats(x, "sum", group=gr)), 
	apply(x, 2L, f(gr, sum)))
expect_equal(
	cbind(row_stats(x, "sum", group=gc)), 
	apply(x, 1L, f(gc, sum)))
expect_equal(
	cbind(col_stats(x, "mean", group=gr)), 
	apply(x, 2L, f(gr, mean)))
expect_equal(
	cbind(row_stats(x, "mean", group=gc)), 
	apply(x, 1L, f(gc, mean)))
expect_equal(
	cbind(col_stats(x, "var", group=gr)), 
	apply(x, 2L, f(gr, var)))
expect_equal(
	cbind(row_stats(x, "var", group=gc)), 
	apply(x, 1L, f(gc, var)))

