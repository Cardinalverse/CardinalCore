require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# peaks find
x1 <- c(0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 2, 2, 3, 0, 1)
x2 <- replace(x1, 11, 3)
x3 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
m1 <- peaks_find(x1)
expect_equal(m1, c(4, 7, 13))
m2 <- peaks_find(x2)
expect_equal(m2, c(4, 7, 11))
m3 <- peaks_find(x3)
expect_equal(m3, 5)

