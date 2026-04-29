require(CardinalCore)
require(testthat)

test_that("qdiff works", {

	expect_equal(qdiff(1:10), diff(1:10))
	expect_equal(qdiff(as.double(1:10)), diff(as.double(1:10)))
	expect_equal(
		qdiff(1:10, units="relative"),
		diff(1:10) / (1:9))
	expect_equal(
		qdiff(as.double(1:10), units="relative"),
		diff(as.double(1:10)) / (1:9))
	expect_equal(
		qdiff(1:10, units="ppm"),
		1e6 * diff(1:10) / (1:9))
	expect_equal(
		qdiff(as.double(1:10), units="ppm"),
		1e6 * diff(as.double(1:10)) / (1:9))
	expect_equal(
		qdiff(as.integer(c(0, NA, NA)), as.integer(c(NA, 0, NA))),
		c(-Inf, Inf, 0))
	expect_equal(
		qdiff(as.double(c(0, NA, NA)), as.double(c(NA, 0, NA))),
		c(-Inf, Inf, 0))
	expect_equal(
		qdiff(c(0, NA, NA), c(NA, 0, NA), units="relative"),
		c(-Inf, Inf, 0))
	expect_equal(
		qdiff(c(0, NA, NA), c(NA, 0, NA), units="ppm"),
		c(-Inf, Inf, 0))
	expect_error(qdiff(LETTERS))

})

test_that("qsort and friends all work", {

	set.seed(1, kind="default")
	u1 <- as.numeric(sample(100L))
	u2 <- as.numeric(sample(101L))
	u3 <- c(0,1,0,1,0,0,3,2,2,2,4,4,8,2,0,0)
	u4 <- c(0,1,NA,1,0,0,3,2,2,NA,4,4,8,2,0,0)

	expect_equal(qselect(u1, 1L), min(u1))
	expect_equal(qselect(u1, 100L), max(u1))
	expect_equal(qselect(u2, 51L), median(u2))
	expect_equal(qselect(u3, 1L), 0)
	expect_equal(qselect(u3, 8L), 1)
	expect_equal(qselect(u3, 9L), 2)
	expect_equal(qselect(u3, 16L), 8)
	expect_equal(qselect(u4, 1L), 0)
	expect_equal(qselect(u4, 8L), 2)
	expect_equal(qselect(u4, 9L), 2)
	expect_equal(qselect(u4, 16L), NA_real_)

	expect_equal(qsort(u1), sort(u1))
	expect_equal(qsort(u2), sort(u2))
	expect_equal(qsort(u3), sort(u3))
	expect_equal(qsort(u4), sort(u4, na.last=TRUE))
	expect_equal(
		qsort(u1, decreasing=TRUE),
		sort(u1, decreasing=TRUE))
	expect_equal(
		qsort(u1, index.return=TRUE),
		sort(u1, index.return=TRUE))
	expect_equal(
		qsort(u1, decreasing=TRUE, index.return=TRUE),
		sort(u1, decreasing=TRUE, index.return=TRUE))

	expect_equal(qmedian(u1), median(u1))
	expect_equal(qmedian(u2), median(u2))
	expect_equal(qmedian(u3), median(u3))
	expect_equal(qmedian(u4, na.rm=TRUE), median(u4, na.rm=TRUE))
	
	expect_equal(qmad(u1), mad(u1))
	expect_equal(qmad(u2), mad(u2))
	expect_equal(qmad(u3), mad(u3))
	expect_equal(qmad(u4, na.rm=TRUE), mad(u3, na.rm=TRUE))

})
