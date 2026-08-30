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
// Summarize peaks in a signal vector

// A vector with peaks (local maxima of k points)
template<Vec U>
struct peaks
{
	U y;
	int k = 5;

	ptrdiff_t ssize() const noexcept { return y.ssize(); }

	auto operator[](ptrdiff_t i) const noexcept
	{
		return is_peak(i) ? y[i] : static_cast<typeof_vec<U>>(0);
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
			if ( diff(y[j], y[i]) >= 0 ) {
				peak = false;
				break;
			}
		}
		for ( ptrdiff_t j = i + 1; j <= i + halfWindow; ++j )
		{
			if ( diff(y[j], y[i]) > 0 ) {
				peak = false;
				break;
			}
		}
		return peak;
	}

	// Get the count of peaks
	ptrdiff_t count() const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			n += is_peak(i);
		return n;
	}

	// Find left endpoint of a peak at i (nearest local minimum)
	ptrdiff_t end_left(ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t end = i > 0 ? i - 1 : 0;
		while ( i > 0 && i >= end - (k / 2) )
		{
			if ( diff(y[i], y[end]) < 0 )
				end = i;
			--i;
		}
		return end;
	}

	// Find right endpoint of a peak at i (nearest local minimum)
	ptrdiff_t end_right(ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t end = i < ssize() - 1 ? i + 1 : ssize() - 1;
		while ( i < ssize() - 1 && i <= end + (k / 2) )
		{
			if ( diff(y[i], y[end]) < 0 )
				end = i;
			++i;
		}
		return end;
	}

	// Find left base of a peak at i (minimum to next higher peak)
	ptrdiff_t base_left(ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t base = i > 0 ? i - 1 : 0;
		for ( ptrdiff_t j = base; j > 0; --j )
		{
			if ( diff(y[j], y[base]) < 0 )
				base = j;
			if ( diff(y[j], y[i]) > 0 )
				break;
		}
		return base;
	}

	// Find right base of a peak at i (minimum to next higher peak)
	ptrdiff_t base_right(ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t base = i < ssize() - 1 ? i + 1 : ssize() - 1;
		for ( ptrdiff_t j = base; j < ssize(); ++j )
		{
			if ( diff(y[j], y[base]) < 0 )
				base = j;
			if ( diff(y[j], y[i]) > 0 )
				break;
		}
		return base;
	}

	// Get indices of peaks and copy into index
	template<Num Index>
	ptrdiff_t index_into(vec<Index> index) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) )
				index[n++] = i;
		}
		return n;
	}

	// Get limits of peaks and copy into output vectors
	template<Num Index>
	ptrdiff_t limits_into(
		vec<Index> index,
		vec<Index> end_left,
		vec<Index> end_right,
		vec<Index> base_left,
		vec<Index> base_right) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				end_left[n] = this->end_left(i);
				end_right[n] = this->end_right(i);
				base_left[n] = this->base_left(i);
				base_right[n] = this->base_right(i);
				++n;
			}
		}
		return n;
	}

};

#endif // CARDINAL_CORE_PEAKS
