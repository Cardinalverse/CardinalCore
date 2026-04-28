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

})
