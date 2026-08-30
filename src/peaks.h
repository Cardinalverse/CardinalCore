#ifndef CARDINAL_CORE_PEAKS
#define CARDINAL_CORE_PEAKS

#include "core.h"
#include "search.h"

//// Utility
//-----------
// Peak utilities

// Get area under the curve by numeric integration
// - Signal y sampled at points x
// - Integral computed within bounds b
// - Estimated using trapezoidal rule
template<Vec U, Vec V>
double trapz(U y, V x, bounds b)
{
	assert(y.ssize() == x.ssize());
	auto y_ = coerce<double>(y);
	auto x_ = coerce<double>(x);
	double sum = 0;
	for ( ptrdiff_t i = b.start; i < b.stop - 1; ++i )
	{
		double dx = x_[i + 1] - x_[i];
		sum += 0.5 * (y_[i + 1] + y_[i]) * dx;
	}
	return sum;
}

//// Vector peaks
//---------------
// Summarize peaks in a signal

// A vector with peaks (local maxima of k points)
template<Vec V>
struct vec_peaks
{
	V x;
	int k = 5;

	ptrdiff_t ssize() const noexcept { return x.ssize(); }

	auto operator[](ptrdiff_t i) const noexcept
	{
		return is_peak(i) ? x[i] : static_cast<typeof_vec<V>>(0);
	}
	
	// Check if element at i is a peak
	// - A peak is > all k/2 points to left
	// - A peak is >= all k/2 points to right
	bool is_peak(ptrdiff_t i) const noexcept
	{
		int halfWindow = k / 2;
		if ( i < halfWindow || i >= ssize() - halfWindow )
			return false;
		bool peak = true;
		for ( ptrdiff_t j = i - 1; j >= i - halfWindow; --j )
		{
			if ( diff(x[j], x[i]) >= 0 ) {
				peak = false;
				break;
			}
		}
		for ( ptrdiff_t j = i + 1; j <= i + halfWindow; ++j )
		{
			if ( diff(x[j], x[i]) > 0 ) {
				peak = false;
				break;
			}
		}
		return peak;
	}

	// Find left endpoint of a peak at i (nearest local minimum)
	ptrdiff_t end_left(ptrdiff_t i) const noexcept
	{
		ptrdiff_t end = i > 0 ? i - 1 : 0;
		while ( i > 0 && i >= end - (k / 2) )
		{
			if ( diff(x[i], x[end]) < 0 )
				end = i;
			--i;
		}
		return end;
	}

	// Find right endpoint of a peak at i (nearest local minimum)
	ptrdiff_t end_right(ptrdiff_t i) const noexcept
	{
		ptrdiff_t end = i < ssize() - 1 ? i + 1 : ssize() - 1;
		while ( i < ssize() - 1 && i <= end + (k / 2) )
		{
			if ( diff(x[i], x[end]) < 0 )
				end = i;
			++i;
		}
		return end;
	}

	// Find left base of a peak at i (minimum to next highest peak)
	ptrdiff_t base_left(ptrdiff_t i) const noexcept
	{
		ptrdiff_t base = i > 0 ? i - 1 : 0;
		for ( ptrdiff_t j = base; j > 0; --j )
		{
			if ( diff(x[j], x[base]) < 0 )
				base = j;
			if ( diff(x[j], x[i]) >= 0 )
				break;
		}
		return base;
	}

	// Find right base of a peak at i (minimum to next highest peak)
	ptrdiff_t base_right(ptrdiff_t i) const noexcept
	{
		ptrdiff_t base = i < ssize() - 1 ? i + 1 : ssize() - 1;
		for ( ptrdiff_t j = base; j < ssize(); ++j )
		{
			if ( diff(x[j], x[base]) < 0 )
				base = j;
			if ( diff(x[j], x[i]) >= 0 )
				break;
		}
		return base;
	}
};

//// Peaks processing
//-------------------
// Process peaks in a signal

// Apply a Callable to each peak in a signal
// - A peak is a local maximum among k points
// - Returns the count of peaks
template<Num Index, UnaryOp Callable, Vec V>
Index peaks_apply(Callable f, const V x, const Index k = 5) noexcept
{
	vec_peaks<V> peaks = {.x = x, .k = k};
	Index count = 0;
	for ( Index i = 0; i < peaks.ssize(); ++i )
	{
		if ( peaks.is_peak(i) )
		{
			f(i);
			++count;
		}
	}
	return count;
}

// Count the number of peaks in a signal
// - A peak is a local maximum among k points
// - Returns the count of peaks
template<Num Index, Vec V>
Index peaks_count(const V x, const Index k = 5) noexcept
{
	return peaks_apply<Index>(sink<Index>{}, x, k);
}

// Find the indices of peaks in a signal
// - Fills index with peak indices (up to index.len)
// - Returns count of peaks (may be > index.len)
template<Num Index, Vec V>
Index peaks_find(vec<Index> index, const V x, const Index k = 5) noexcept
{
	return peaks_apply<Index>(sink<Index>{index}, x, k);
}

#endif // CARDINAL_CORE_PEAKS
