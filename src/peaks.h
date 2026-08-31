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

	// Sum of peak at i
	template<Num T = double>
	T sum(ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
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

	// Area of peak at i
	// - Integrates signal y sampled at x
	// - Computed using trapezoidal rule
	template<Num T = double, Vec V>
	T area(V x, ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		return area(left_end(i), right_end(i));
	}

	// Area of signal on closed interval [lo, hi]
	// - Integrates signal y sampled at x
	// - Computed using trapezoidal rule
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

	// Prominence of peak at i (height above lowest contour line)
	// - Lowest contour line is the higher of its bases
	// - Window wlen gives number of points to search
	// - Window is centered on the peak at i
	template<Num T = double>
	T prominence(ptrdiff_t i, double wlen = 0) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t wl = wlen > 0 ? wlen : ssize();
		return prominence(i, left_base(i, wl), right_base(i, wl));
	}

	// Prominence of peak at i (given bases at lo and hi)
	template<Num T = double>
	T prominence(ptrdiff_t i, ptrdiff_t lo, ptrdiff_t hi) const noexcept
	{
		T ylo = coerce_cast<T>(y[lo]);
		T yhi = coerce_cast<T>(y[hi]);
		T ymax = ylo;
		for ( ptrdiff_t i = lo; i <= hi; ++i )
			ymax = y[i] > ymax ? y[i] : ymax;
		return y[i] - (ylo > yhi ? ylo : yhi);
	}

	// Width of peak at i (at a fraction of max)
	template<Num T = double, Vec V>
	T width(V x, ptrdiff_t i, double fmax = 0.5) const noexcept
	{
		assert(is_peak(i));
		return right_ips<T>(x, i, fmax) - left_ips<T>(x, i, fmax);
	}

	// Centroid of peak at i
	template<Num T = double, Vec V>
	T centroid(V x, ptrdiff_t i) const noexcept
	{
		assert(is_peak(i));
		return centroid(x, left_end(i), right_end(i));
	}

	// Centroid of peak at i
	template<Num T = double, Vec V>
	T centroid(V x, ptrdiff_t lo, ptrdiff_t hi) const noexcept
	{
		T sxy = 0;
		T sy = 0;
		for ( ptrdiff_t i = lo; i <= hi; ++i )
		{
			sxy += y[i] * x[i];
			sy += y[i];
		}
		return sxy / sy;
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

	// Find left base of a peak at i (minimum to next higher peak)
	// - Window wlen gives number of points to search
	// - Window is centered on the peak at i
	ptrdiff_t left_base(ptrdiff_t i, ptrdiff_t wlen = 0) const noexcept
	{
		assert(is_peak(i));
		ptrdiff_t wlo = wlen > 0 ? i - (wlen / 2) : 0;
		ptrdiff_t base = i > 0 ? i - 1 : 0;
		for ( ptrdiff_t j = base; j > 0 && j >= wlo; --j )
		{
			if ( y[j] < y[base] )
				base = j;
			if ( y[j] > y[i] )
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
		ptrdiff_t whi = wlen > 0 ? i + (wlen / 2) : ssize() - 1;
		ptrdiff_t base = i < ssize() - 1 ? i + 1 : ssize() - 1;
		for ( ptrdiff_t j = base; j < ssize() - 1 && j <= whi; ++j )
		{
			if ( y[j] < y[base] )
				base = j;
			if ( y[j] > y[i] )
				break;
		}
		return base;
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
		ptrdiff_t wl = wlen > 0 ? wlen : ssize();
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				left_base[n] = this->left_base(i, wl);
				right_base[n] = this->right_base(i, wl);
				prominences[n] = prominence<T>(i, left_base[n], right_base[n]);
				++n;
			}
		}
		return n;
	}

	// Get areas of peaks and copy into output vectors
	template<Vec V, Num Index, Num T = double>
	ptrdiff_t areas_into(
		vec<Index> index,
		vec<Index> left_end,
		vec<Index> right_end,
		vec<T> areas,
		const V x) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				left_end[n] = this->left_end(i);
				right_end[n] = this->right_end(i);
				areas[n] = area<T>(x, left_end[n], right_end[n]);
				++n;
			}
		}
		return n;
	}

	// Get widths of peaks and copy into output vectors
	template<Vec V, Num Index, Num T = double>
	ptrdiff_t widths_into(
		vec<Index> index,
		vec<T> left_ips,
		vec<T> right_ips,
		vec<T> widths,
		const V x,
		const double fmax = 0.5) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				left_ips[n] = this->left_ips<T>(x, i, fmax);
				right_ips[n] = this->right_ips<T>(x, i, fmax);
				widths[n] = right_ips[n] - left_ips[n];
				++n;
			}
		}
		return n;
	}

	// Get centroids of peaks and copy into output vectors
	template<Vec V, Num Index, Num T = double>
	ptrdiff_t centroids_into(
		vec<Index> index,
		vec<Index> left_end,
		vec<Index> right_end,
		vec<T> centroids,
		const V x) const noexcept
	{
		ptrdiff_t n = 0;
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
		{
			if ( is_peak(i) ) {
				index[n] = i;
				left_end[n] = this->left_end(i);
				right_end[n] = this->right_end(i);
				centroids[n] = centroid<T>(x, left_end[n], right_end[n]);
				++n;
			}
		}
		return n;
	}

};

#endif // CARDINAL_CORE_PEAKS
