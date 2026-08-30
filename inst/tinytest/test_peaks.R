require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# helper
named <- function(x) setNames(x, seq_along(x))

# peaks find
x1 <- c(0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 2, 2, 3, 0, 1)
x2 <- replace(x1, 11, 3)
x3 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
m1 <- peaks_find(x1)
m2 <- peaks_find(x2)
m3 <- peaks_find(x3)
expect_equal(m1, c(4, 7, 13))
expect_equal(m2, c(4, 7, 11))
expect_equal(m3, 5)

b1 <- peaks_find_limits(x1)
b2 <- peaks_find_limits(x2)
b3 <- peaks_find_limits(x3)

x4 <- c(rep(0, 5), rep(1, 5), 2, rep(0, 5))
peaks_find_limits(x4)

