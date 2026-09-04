require(CardinalCore, quietly=TRUE)
require(tinytest, quietly=TRUE)

# filt1
set.seed(1)
y1 <- runif(10)
expect_equal(
	filt1_mean(y1, k=5L),
	roll_apply(y1, k=5L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_mean(y1, k=7L),
	roll_apply(y1, k=7L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_conv(y1, w=rep(1, 5)),
	roll_apply(y1, k=5L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_conv(y1, w=rep(1, 7)),
	roll_apply(y1, k=7L, FUN=mean, na.rm=TRUE))

# filt1 with NAs
y2 <- replace(y1, 5L, NA)
expect_equal(
	filt1_mean(y2, k=5L),
	roll_apply(y2, k=5L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_mean(y2, k=7L),
	roll_apply(y2, k=7L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_conv(y2, w=rep(1, 5)),
	roll_apply(y2, k=5L, FUN=mean, na.rm=TRUE))
expect_equal(
	filt1_conv(y2, w=rep(1, 7)),
	roll_apply(y2, k=7L, FUN=mean, na.rm=TRUE))

