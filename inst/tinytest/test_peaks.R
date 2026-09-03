require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# helper
named <- function(x) setNames(x, seq_along(x))

# peaks find
y1 <- c(0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 2, 2, 3, 0, 1)
y2 <- replace(y1, 11, 3)
y3 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
i1 <- peaks_find(y1)
i2 <- peaks_find(y2)
i3 <- peaks_find(y3)
expect_equal(i1, c(4, 7, 13))
expect_equal(i2, c(4, 7, 11))
expect_equal(i3, 5)

p1 <- peaks_prominences(y1)
p2 <- peaks_prominences(y2)
p3 <- peaks_prominences(y3)
expect_equal(p1$left_base, c(3, 6, 10))
expect_equal(p1$right_base, c(5, 10, 14))
expect_equal(p1$prominence, c(1, 1, 3))
expect_equal(p2$prominence, c(1, 1, 3))
expect_equal(p3$prominence, 1)

peaks_summary(y1)
peaks_summary(y2)
peaks_summary(y3)

