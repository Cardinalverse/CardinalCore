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

# test

path <- "/Volumes/Local/Data/public/pride/PXD001283/HR2MSI mouse urinary bladder S096.imzML"
mzml <- CardinalIO::parseImzML(path, ibd=TRUE)

i <- 100
y <- mzml$ibd$intensity[[i]]
x <- mzml$ibd$mz[[i]]
k <- 5
n <- 1001
noise <- roll_apply(y - filt1_mean(y, k=k), k=n, FUN=mad)
noise_list <- roll(y - filt1_mean(y, k=k), k=n)
snr <- y/noise
threshold <- 12
par(mfrow=c(2,1))
plot(y ~ x, type="l")
points(y[snr > threshold] ~ x[snr > threshold], col="red")
plot(snr ~ x, type="h")
abline(h=6, col="blue")
cbind(y,snr,noise)

