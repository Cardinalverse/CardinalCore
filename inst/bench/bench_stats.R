require(CardinalCore, quietly=TRUE)

# vec_stats
set.seed(1)
n <- 1e7
x <- runif(n)
g <- sample(6L, n, replace=TRUE)
f <- factor(g)
size <- 1e3
m <- rep.int(size, n %/% size)

bench::mark(
	vec_stats(x, "sum", lengths=m, num.threads=1),
	vec_stats(x, "sum", lengths=m, num.threads=2),
	vec_stats(x, "sum", lengths=m, num.threads=4))

bench::mark(
	vec_stats(x, "sum", g, lengths=m, num.threads=1),
	vec_stats(x, "sum", g, lengths=m, num.threads=2),
	vec_stats(x, "sum", g, lengths=m, num.threads=4))

bench::mark(
	vec_stats(x, "sum", f, lengths=m, num.threads=1),
	vec_stats(x, "sum", f, lengths=m, num.threads=2),
	vec_stats(x, "sum", f, lengths=m, num.threads=4))

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

