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

# test case 1

path <- "/Volumes/Local/Data/public/pride/PXD001283/HR2MSI mouse urinary bladder S096.imzML"
mzml <- CardinalIO::parseImzML(path, ibd=TRUE)

library(S4Vectors)
library(IRanges)
library(tinyplot)
tinytheme("dark")
spec <- function(i) {
	cbind(mz=mzml$ibd$mz[[i]], intensity=mzml$ibd$intensity[[i]])
}
s <- spec(1)
sdf <- as.data.frame(s)
sDF <- DataFrame(s)
peakPick(s)
peakPick(sdf)
peakPick(sDF)
mzr <- range(sdf$mz)
mzs <- seq(from=mzr[1], to=mzr[2], by=0.1)

bench::mark(peakPick(s))
bench::mark(peakPick(sdf))
bench::mark(peakPick(sDF))

group_by_ref(sdf$mz, ref=mzs, tolerance=0.05)
group_by_ref(sdf, ref=list(mz=mzs), tolerance=0.05)

s |> peakPick() |> peakAlign(ref=list(mz=mzs), tolerance=0.05)
sdf |> peakPick() |> peakAlign(ref=list(mz=mzs), tolerance=0.05)
sDF |> peakPick() |> peakAlign(ref=list(mz=mzs), tolerance=0.05)

bench::mark(s |> peakPick() |> peakAlign(ref=list(centroid=mzs), ppm=80))
bench::mark(sdf |> peakPick() |> peakAlign(ref=list(centroid=mzs), ppm=80))
bench::mark(sDF |> peakPick() |> peakAlign(ref=list(centroid=mzs), ppm=80))

s |> peakPick() |> peakAlign(ref=list(mz=mzs), tolerance=0.045)
s |> peakPick() |> peakAlign(ref=list(centroid=mzs), tolerance=0.045)

intensity <- function(i) mzml$ibd$intensity[[i]]
mz <- function(i) mzml$ibd$mz[[i]]
process <- function(i)
{
	if ( i %% 1000L == 0L ) message(i, "/", length(mzml$ibd$mz))
	as.data.frame(peaks_summary(intensity(i), mz(i)))
}
head(process(505), n=20)

intensities <- as.list(mzml$ibd$intensity)
mzs <- as.list(mzml$ibd$mz)
process <- function(i)
{
	if ( i %% 1000L == 0L ) message(i, "/", length(mzs))
	as.data.frame(peaks_summary(intensities[[i]], mzs[[i]]))
}
head(process(505), n=20)

system.time(peaks <- lapply(seq_along(mzml$ibd$mz), process))


head(p$max / matter::estnoise_diff(y)[1L])

# test case 2

path <- "/Volumes/Local/Data/private/scratch/timsdata/input.imzML"
mzml <- CardinalIO::parseImzML(path, ibd=TRUE, extraArrays=c(mobility="MS:1003006"))

spec <- function(i) {
	data.frame(
		mz=mzml$ibd$mz[[i]],
		mobility=mzml$ibd$extra$mobility[[i]],
		intensity=mzml$ibd$intensity[[i]])
}
str(mzml$run$spectrumList)

i <- 9000
sdf <- spec(i)
sDF <- DataFrame(sdf)
mzr <- range(sdf$mz)
mobr <- range(sdf$mobility)
mzs <- seq(from=floor(mzr[1]), to=ceiling(mzr[2]), by=0.05)
mobs <- seq(from=floor(mobr[1]), to=ceiling(mobr[2]), length.out=10)
dmz <- max(diff(mzs)) / 2
dmobs <- max(diff(mobs)) / 2
sDF2 <- peakAlign(sDF,
	ref=list(mobility=mobs, mz=mzs),
	tolerance=c(mobility=dmobs, mz=dmz))
print(sDF2)

hits <- msearch(s[c("mz", "mobility")], list(mzs, mobs), tolerance=c(dmz, dmobs))

xmzs <- peaks_agg(s$intensity, s$mz, mzs, stat="sum", tolerance=dmz)
xmobs <- peaks_agg(s$intensity, s$mobility, mobs, stat="sum", tolerance=dmobs)

require(tinyplot)
tinytheme("dark")
xlim <- c(800, 810)
xlim <- c(806.55, 806.59)
plt(s$mz, s$mobility, cex=0.1)
plt(s$mz, s$mobility, cex=s$intensity, xlim=xlim)
plt(s$mz, s$mobility, cex=0.1, alpha=0.5, by=log1p(s$intensity), fill="by", pch=21, palette=hcl.colors(100, palette="viridis", rev=TRUE))
plt(s$mz, s$mobility, cex=0.5, by=log1p(s$intensity), fill="by", pch=21, xlim=xlim, ylim=c(1.35, 1.5), palette=hcl.colors(100, palette="viridis", rev=TRUE))
plt(s$mz, s$intensity, by=s$mobility, cex=0.5, xlim=xlim)
plt(mzs, xmzs, type='l', xlim=xlim)
plt(mobs, xmobs, type='l')
s[23400:23459,]

gmz <- bsearch(s$mz, mzs, tolerance=dmz)
gmob <- bsearch(s$mobility, mobs, tolerance=dmobs)
g <- cbind(gmz, gmob)
m <- kdsearch_agg(g, g, seq_len(nrow(g)), stat="min")
xg <- group_stats(stream_stats(s$intensity, stat="mean"), m)
mzg <- group_stats(stream_stats(s$mz, stat="mean"), m)
mobg <- group_stats(stream_stats(s$mobility, stat="mean"), m)
sp <- data.frame(mz=mzg, mobility=mobg, intensity=xg, n=attr(xg, "nobs"))

