require(CardinalCore, quietly=TRUE)

# vec_stats
set.seed(1)
n1 <- 9000L
n2 <- 1000L
xs <- replicate(n1, runif(n2), simplify=FALSE)
gs <- replicate(n1, sample(6L, n2, replace=TRUE), simplify=FALSE)
l <- lengths(xs)
x <- unlist(xs)
g <- unlist(gs)
f <- factor(g)

bench::mark(
	vec_stats(x, "mean", num.threads=0),
	vec_stats(x, "mean", num.threads=1),
	vec_stats(x, "mean", num.threads=2),
	vec_stats(x, "mean", num.threads=4),
	vec_stats(x, "mean", num.threads=8))

bench::mark(
	vec_stats(x, "mean", f, num.threads=0),
	vec_stats(x, "mean", f, num.threads=1),
	vec_stats(x, "mean", f, num.threads=2),
	vec_stats(x, "mean", f, num.threads=4),
	vec_stats(x, "mean", f, num.threads=8))

bench::mark(
	vec_stats(x, "mean", f, sizes=l, collapse=TRUE, num.threads=0),
	vec_stats(x, "mean", f, sizes=l, collapse=TRUE, num.threads=1),
	vec_stats(x, "mean", f, sizes=l, collapse=TRUE, num.threads=2),
	vec_stats(x, "mean", f, sizes=l, collapse=TRUE, num.threads=4),
	vec_stats(x, "mean", f, sizes=l, collapse=TRUE, num.threads=8))

# col_stats and row_stats
set.seed(1)
nr <- 7000L
nc <- 5000L
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
x[sample(100)] <- NA_real_
gr <- sample(6L, nr, replace=TRUE)
gc <- sample(6L, nc, replace=TRUE)
fr <- factor(gr)
fc <- factor(gc)
tx <- t(x)

bench::mark(
	col_stats(x, "mean", num.threads=0),
	col_stats(x, "mean", num.threads=1),
	col_stats(x, "mean", num.threads=2),
	col_stats(x, "mean", num.threads=4),
	col_stats(x, "mean", num.threads=8))

bench::mark(
	row_stats(x, "mean", num.threads=0),
	row_stats(x, "mean", num.threads=1),
	row_stats(x, "mean", num.threads=2),
	row_stats(x, "mean", num.threads=4),
	row_stats(x, "mean", num.threads=8))

bench::mark(
	col_stats(x, "mean", fr, num.threads=0),
	col_stats(x, "mean", fr, num.threads=1),
	col_stats(x, "mean", fr, num.threads=2),
	col_stats(x, "mean", fr, num.threads=4),
	col_stats(x, "mean", fr, num.threads=8))

bench::mark(
	row_stats(x, "mean", fc, num.threads=0),
	row_stats(x, "mean", fc, num.threads=1),
	row_stats(x, "mean", fc, num.threads=2),
	row_stats(x, "mean", fc, num.threads=4),
	row_stats(x, "mean", fc, num.threads=8))

