require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# setup
set.seed(1, kind="default")
u1 <- as.numeric(sample(100L))
u2 <- as.numeric(sample(101L))
u3 <- c(0,1,0,1,0,0,3,2,2,2,4,4,8,2,0,0)
u4 <- c(0,1,NA,1,0,0,3,2,2,NA,4,4,8,2,0,0)

# qselect
expect_equal(qselect(u1, c(1L, 100L)), c(min(u1), max(u1)))
expect_equal(qselect(u2, c(1L, 51L)), c(min(u1), median(u2)))
expect_equal(qselect(u3, c(1L, 8L, 9L, 16L)), c(0, 1, 2, 8))
expect_equal(qselect(u4, c(1L, 8L, 9L, 16L)), c(0, 2, 2, NA_real_))
expect_error(qselect(LETTERS, 1L))

# qorder
expect_equal(qorder(u1), order(u1))
expect_equal(qorder(u2), order(u2))
expect_equal(qorder(u3), order(u3))
expect_equal(qorder(u4), order(u4))
expect_error(qorder(LETTERS))

# qmedian
expect_equal(qmedian(u1), median(u1))
expect_equal(qmedian(u2), median(u2))
expect_equal(qmedian(u3), median(u3))
expect_equal(qmedian(u4), median(u4, na.rm=TRUE))
expect_error(qmedian(LETTERS))

# qmad
expect_equal(qmad(u1), mad(u1))
expect_equal(qmad(u2), mad(u2))
expect_equal(qmad(u3), mad(u3))
expect_equal(qmad(u4), mad(u3, na.rm=TRUE))
expect_error(qmad(LETTERS))

