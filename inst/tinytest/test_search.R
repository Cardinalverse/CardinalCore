require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# bsearch on integers
q <- c(1L, 4L, 6L, 99L, 102L)
x <- c(1L, 2L, 3L, 4L, 5L, 8L, 9L, 100L, 101L, 102L)
expect_equal(bsearch(q, x), c(1, 4, NA, NA, 10))
expect_equal(bsearch(q, x, tolerance=Inf), c(1, 4, 5, 8, 10))
expect_equal(bsearch((-1):1, integer(), tolerance=0.1), rep(NA_integer_, 3L))

# bsearch on doubles
q <- c(1.11, 3.0, 3.33, 5.0, 5.1)
x <- c(1.11, 2.22, 3.33, 4.0, 5.0)
expect_equal(bsearch(q, x), c(1, NA, 3, 5, NA))
expect_equal(bsearch(q, x, tolerance=Inf), c(1, 3, 3, 5, 5))
expect_equal(bsearch(3.0, x, tolerance=0), NA_integer_)
expect_equal(bsearch(3.0, x, tolerance=0.5), 3)
expect_equal(bsearch(3.0, x, tolerance=Inf), 3)
expect_equal(bsearch(3.0, x, tolerance=0.11, relative_to="query"), NA_integer_)
expect_equal(bsearch(3.0, x, tolerance=0.11, relative_to="table"), 3L)
expect_equal(bsearch(3.0, x, tolerance=Inf, relative_to="query"), 3)
expect_equal(bsearch(3.0, x, tolerance=Inf, relative_to="table"), 3)
expect_equal(bsearch(3.0, x, tolerance=0.111, relative_to="query"), 3)
expect_equal(bsearch(3.0, x, tolerance=0.111, relative_to="table"), 3)
expect_equal(bsearch((-1):1, numeric(), tolerance=0.1), rep(NA_integer_, 3L))

# kdtree build
d1 <- data.frame(
	x=c(2,5,9,4,8,7,9,8,9,6,3,1,9,2,8),
	y=c(3,4,6,7,1,2,4,4,7,3,4,6,5,1,7),
	z=c(3,2,7,9,5,6,1,2,8,1,5,8,3,3,6))
d2 <- data.frame(
	x=c(1,2,3,4,1,2,3,5,9,3,1),
	y=c(1,2,3,4,1,2,3,4,7,3,1))
i1 <- seq_len(nrow(d1))
i2 <- seq_len(nrow(d2))
t1 <- kdtree(d1)
t2 <- kdtree(d2)
ns1 <- cbind(t1$left + 1L, t1$right + 1L)
expect_equal(t1$root + 1L, 6L)
expect_equal(t1$root + 1L, which(!i1 %in% sort(ns1, na.last=NA)))

# kdtree range search
ks1a <- kdsearch(c(2,3,3), d1, tolerance=2)
ks1b <- kdsearch(c(7,2,6), d1, tolerance=2)
ks1c <- kdsearch(rbind(c(2,3,3), c(7,2,6)), d1, tolerance=2)
expect_equal(ks1a$index, c(1L, 11L, 14L))
expect_equal(ks1b$index, c(5L, 6L))
expect_equal(ks1c$index, c(ks1a$index, ks1b$index))
ks2a <- kdsearch(c(2,3,3), d1, tolerance=c(2,2,4))
ks2b <- kdsearch(c(7,2,6), d1, tolerance=c(2,2,4))
ks2c <- kdsearch(rbind(c(2,3,3), c(7,2,6)), d1, tolerance=c(2,2,4))
expect_equal(ks2a$index, c(1L, 11L, 14L))
expect_equal(ks2b$index, c(2L, 5L, 6L, 8L))
expect_equal(ks2c$index, c(ks2a$index, ks2b$index))

# kdtree knn search
knnsearch(c(2,3,3), d1, k=3)
knnsearch(c(7,2,6), d1, k=3)
knnsearch(rbind(c(2,3,3), c(7,2,6)), d1, k=3)

