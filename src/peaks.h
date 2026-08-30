#ifndef CARDINAL_CORE_PEAKS
#define CARDINAL_CORE_PEAKS

#include "core.h"
#include "search.h"

//// Peaks vector
//---------------
// Summarize peaks in a signal

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
		int hk = k / 2;
		if ( i < hk || i >= ssize() - hk )
			return false;
		bool peak = true;
		for ( ptrdiff_t j = i - 1; j >= i - hk; --j )
		{
			if ( diff(y[j], y[i]) >= 0 ) {
				peak = false;
				break;
			}
		}
		for ( ptrdiff_t j = i + 1; j <= i + hk; ++j )
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
	ptrdiff_t left_end(ptrdiff_t i) const noexcept
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
	ptrdiff_t right_end(ptrdiff_t i) const noexcept
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

	// Sum of peak at i
	template<Num T = double>
	T sum(ptrdiff_t i) const noexcept {
		return sum(left_end(i), right_end(i));
	}

	// Sum of signal on closed interval [lo, hi]
	template<Num T = double>
	T sum(ptrdiff_t lo, ptrdiff_t hi) const noexcept
	{
		T sum = 0;
		for ( ptrdiff_t i = lo; i <= hi; ++i )
			sum += y[i];
		return sum;
	}

	// Area of peak at i (trapezoidal integration)
	template<Num T = double, Vec V>
	T area(V x, ptrdiff_t i) const noexcept {
		return area(left_end(i), right_end(i));
	}

	// Area of signal on closed interval [lo, hi] (trapezoidal integration)
	template<Num T = double, Vec V>
	T area(V x, ptrdiff_t lo, ptrdiff_t hi) const noexcept
	{
		assert(x.ssize() == y.ssize());
		T area = 0;
		auto y_ = coerce<T>(y);
		auto x_ = coerce<T>(x);
		for ( ptrdiff_t i = lo; i < hi; ++i )
		{
			T dx = x_[i + 1] - x_[i];
			area += 0.5 * (y_[i + 1] + y_[i]) * dx;
		}
		return area;
	}

	// Find left base of a peak at i (minimum to next higher peak)
	// - Window wlen gives number of points to search
	// - Window is centered on the peak at i
	ptrdiff_t left_base(ptrdiff_t i, ptrdiff_t wlen = 0) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t w = wlen > 0 ? wlen : ssize();
		ptrdiff_t wmin = i - (w / 2);
		ptrdiff_t base = i > 0 ? i - 1 : 0;
		for ( ptrdiff_t j = base; j > 0 && j > wmin; --j )
		{
			if ( diff(y[j], y[base]) < 0 )
				base = j;
			if ( diff(y[j], y[i]) > 0 )
				break;
		}
		return base;
	}

	// Find right base of a peak at i (minimum to next higher peak)
	// - Window wlen gives number of points to search
	// - Window is centered on the peak at i
	ptrdiff_t right_base(ptrdiff_t i, ptrdiff_t wlen = 0) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t w = wlen > 0 ? wlen : ssize();
		ptrdiff_t wmax = i + (w / 2);
		ptrdiff_t base = i < ssize() - 1 ? i + 1 : ssize() - 1;
		for ( ptrdiff_t j = base; j < ssize() - 1 && j < wmax; ++j )
		{
			if ( diff(y[j], y[base]) < 0 )
				base = j;
			if ( diff(y[j], y[i]) > 0 )
				break;
		}
		return base;
	}

	// Width of peak at i (at a fraction of max)
	template<Num T = double>
	T prominence(ptrdiff_t i, double wlen = 0) const noexcept
	{
		ptrdiff_t w = wlen > 0 ? wlen : ssize();
		return prominence(i, left_base(i, w), right_base(i, w));
	}

	// Width of peak at i given bases at [lo, hi]
	template<Num T = double>
	T prominence(ptrdiff_t i, ptrdiff_t lo, ptrdiff_t hi) const noexcept
	{
		T ylo = coerce_cast<T>(y[lo]);
		T yhi = coerce_cast<T>(y[hi]);
		return y[i] - (ylo > yhi ? ylo : yhi);
	}

	// Find left intersection at a fraction of the max of peak i
	template<Num T = double, Vec V>
	T left_ips(V x, ptrdiff_t i, double fmax = 0.5) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t j = i > 0 ? i - 1 : 0;
		while ( j >= 0 )
		{
			T fheight = fmax * y[i];
			if ( y[j] < fheight )
			{
				T t = (fheight - y[j]) / (y[j + 1] - y[j]);
				return x[j] + (t * (x[j + 1] - x[j]));
			}
			else
				--j;
		}
		return x[j];
	}

	// Find right intersection at a fraction of the max of peak i
	template<Num T = double, Vec V>
	T right_ips(V x, ptrdiff_t i, double fmax = 0.5) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t j = i < ssize() - 1 ? i + 1 : ssize() - 1;
		while ( j < ssize() )
		{
			T fheight = fmax * y[i];
			if ( y[j] < fheight )
			{
				T t = (fheight - y[j - 1]) / (y[j] - y[j - 1]);
				return x[j - 1] + (t * (x[j] - x[j - 1]));
			}
			else
				++j;
		}
		return x[j];
	}

	// Width of peak at i (at a fraction of max)
	template<Num T = double, Vec V>
	T width(V x, ptrdiff_t i, double fmax = 0.5) const noexcept
	{
		return right_ips<T>(x, i, fmax) - left_ips<T>(x, i, fmax);
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

	// Get sums of peaks and copy into output vectors
	template<Num Index, Num T = double>
	ptrdiff_t sums_into(
		vec<Index> index,
		vec<Index> left_end,
		vec<Index> right_end,
		vec<T> sums) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) )
			{
				index[n] = i;
				left_end[n] = this->left_end(i);
				right_end[n] = this->right_end(i);
				sums[n] = sum<T>(left_end[n], right_end[n]);
				++n;
			}
		}
		return n;
	}

	// Get prominences of peaks and copy into output vectors
	template<Num Index, Num T = double>
	ptrdiff_t prominences_into(
		vec<Index> index,
		vec<Index> left_base,
		vec<Index> right_base,
		vec<T> prominences,
		const ptrdiff_t wlen = 0) const noexcept
	{
		ptrdiff_t n = 0;
		ptrdiff_t w = wlen > 0 ? wlen : ssize();
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				left_base[n] = this->left_base(i, w);
				right_base[n] = this->right_base(i, w);
				prominences[n] = prominence<T>(i, left_base[n], right_base[n]);
				++n;
			}
		}
		return n;
	}

};

#endif // CARDINAL_CORE_PEAKS
