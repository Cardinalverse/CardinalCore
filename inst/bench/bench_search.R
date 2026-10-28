require(CardinalCore, quietly=TRUE)

# bsearch on integers
q <- runif(1e4)
x <- sort(runif(1e6))

bench::mark(
	bsearch(q, x, tolerance=0.1, num.threads=1),
	bsearch(q, x, tolerance=0.1, num.threads=2),
	bsearch(q, x, tolerance=0.1, num.threads=4),
	bsearch(q, x, tolerance=0.1, num.threads=8))

bench::mark(
	bsearch_agg(q, x, x, tolerance=0.1, num.threads=1),
	bsearch_agg(q, x, x, tolerance=0.1, num.threads=2),
	bsearch_agg(q, x, x, tolerance=0.1, num.threads=4),
	bsearch_agg(q, x, x, tolerance=0.1, num.threads=8))


