require(CardinalCore, quietly=TRUE)

# vec_stats
set.seed(1)
n1 <- 1e4
n2 <- 1000L
xs <- replicate(n1, runif(n2), simplify=FALSE)
gs <- replicate(n1, sample(6L, n2, replace=TRUE), simplify=FALSE)
l <- lengths(xs)
x <- unlist(xs)
g <- unlist(gs)
f <- factor(g)

bench::mark(
	vec_stats(x, "sum", l, num.threads=1),
	vec_stats(x, "sum", l, num.threads=2),
	vec_stats(x, "sum", l, num.threads=4))

bench::mark(
	vec_stats(x, "sum", l, g, num.threads=1),
	vec_stats(x, "sum", l, g, num.threads=2),
	vec_stats(x, "sum", l, g, num.threads=4))

bench::mark(
	vec_stats(x, "sum", l, f, num.threads=1),
	vec_stats(x, "sum", l, f, num.threads=2),
	vec_stats(x, "sum", l, f, num.threads=4))

# col_stats and row_stats
set.seed(1)
nr <- 1e6
nc <- 32L
x <- matrix(runif(nr * nc), nrow=nr, ncol=nc)
x[sample(100)] <- NA_real_
gr <- sample(6L, nr, replace=TRUE)
gc <- sample(6L, nc, replace=TRUE)
fr <- factor(gr)
fc <- factor(gc)

bench::mark(
	col_stats(x, "sum", num.threads=1),
	col_stats(x, "sum", num.threads=2),
	col_stats(x, "sum", num.threads=4))

bench::mark(
	col_stats(x, "sum", gr, num.threads=1),
	col_stats(x, "sum", gr, num.threads=2),
	col_stats(x, "sum", gr, num.threads=4))

bench::mark(
	col_stats(x, "sum", fr, num.threads=1),
	col_stats(x, "sum", fr, num.threads=2),
	col_stats(x, "sum", fr, num.threads=4))

