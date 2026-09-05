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

# test

path <- "/Volumes/Local/Data/public/pride/PXD001283/HR2MSI mouse urinary bladder S096.imzML"
mzml <- CardinalIO::parseImzML(path, ibd=TRUE)

i <- 2000
y <- mzml$ibd$intensity[[i]]
x <- mzml$ibd$mz[[i]]
bench::mark(p <- peaks_summary(y, x))
p <- as.data.frame(p)
head(p, n=20)

intensity <- function(i) mzml$ibd$intensity[[i]]
mz <- function(i) mzml$ibd$mz[[i]]
process <- function(i)
{
	if ( i %% 1000L == 0L ) message(i, "/", length(mzml$ibd$mz))
	as.data.frame(peaks_summary(intensity(i), mz(i)))
}
head(process(505), n=20)

system.time(peaks <- lapply(seq_along(mzml$ibd$mz), process))

head(p$max / matter::estnoise_diff(y)[1L])

# test 2

path <- "/Volumes/Local/Data/private/scratch/timsdata/output.d/input.imzML"
mzml <- CardinalIO::parseImzML(path, ibd=TRUE, extraArrays=c(mobility="MS:1003006"))

mzs <- readBin(memDecompress(mzml$ibd$mz[[1]]),
	what="double", n=40309)
mobs <- readBin(memDecompress(mzml$ibd$extra$mobility[[1]]),
	what="double", n=40309)
ints <- readBin(memDecompress(mzml$ibd$intensity[[1]]),
	what="double", n=40309)
df <- data.frame(mz=mzs, mobility=mobs, intensity=ints)
plot(mobility ~ mz, data=df, cex=df$intensity / max(df$intensity))

