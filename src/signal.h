#ifndef CARDINAL_CORE_SIGNAL
#define CARDINAL_CORE_SIGNAL

#include "core.h"
#include "search.h"

//// Utility
//----------
// Signal processing utilities

template<Num L, Num R>
L min2(const L x, const R y) noexcept { return x < y ? x : y; }

template<Num L, Num R>
L max2(const L x, const R y) noexcept { return x > y ? x : y; }

// Clamp x to [lo, hi]
template<Num T, Num L, Num R>
T clamp(const T x, const L lo, const R hi) noexcept
{
	return min2(max2(x, lo), hi);
}

//// Shifted signal
//-----------------
// Vector with lags or leads

template<Vec U, Num T = double>
struct vec_shifted
{
	U y;
	int shift;

	ptrdiff_t ssize() const noexcept { return y.ssize(); }

	T operator[](const ptrdiff_t i) const noexcept
	{
		ptrdiff_t j = i + shift;
		if ( 0 <= j && j < ssize() )
			return y[j];
		else
			return na_value<T>();
	}
};

template<Num T = double, Vec U>
auto shift(const U y, int shift) noexcept
{
	return vec_shifted<U,T>{y, shift};
}

//// Convolved signal
//--------------------
// Vector convolved with a window

template<Vec U, Vec W, Num T = double>
struct vec_convolved
{
	U y;
	W w;

	ptrdiff_t ssize() const noexcept { return y.ssize(); }

	T operator[](const ptrdiff_t i) const noexcept
	{
		T sw = 0;
		T swy = 0;
		ptrdiff_t iw = 0;
		ptrdiff_t iy = i - (w.ssize() / 2);
		while ( iw < w.ssize() )
		{
			ptrdiff_t j = clamp(iy, 0, y.ssize() - 1);
			if ( !is_na(w[iw]) && !is_na(y[j]) )
			{
				sw += w[iw];
				swy += w[iw] * y[j];
			}
			++iw;
			++iy;
		}
		return swy / sw;
	}
};

template<Num T = double, Vec U, Vec W>
auto convolve(const U y, const W w) noexcept
{
	return vec_convolved<U,W,T>{y, w};
}

//// Kernel vectors
//------------------
// Windows for convolution filters

template<Num T = double>
struct gaussian
{
	ptrdiff_t len;
	T sigma;

	ptrdiff_t ssize() const noexcept { return len; }

	T operator[](ptrdiff_t i) const noexcept
	{
		ptrdiff_t x = i - (len / 2);
		return std::exp(-(x * x) / (2 * sigma * sigma));
	}
};

//// 1D Filters
//-------------
// Filter 1D signals

// Mean filter
template<Num T, Vec V>
void filt1_mean(vec<T> dst, V src, int k = 5) noexcept
{
	assert(dst.len == src.ssize());
	int hk = k / 2;
	int n = 0;
	T sum = 0;
	// initialize sliding window
	for ( ptrdiff_t i = 0; i <= hk && i < dst.len; ++i )
	{
		if ( is_na(src[i]) )
			continue;
		// initialize first element and out-of-signal elements
		if ( i == 0 )
		{
			sum += (hk + 1) * src[0];
			n += (hk + 1);
		}
		// initialize remaining elements
		else
		{
			sum += src[i];
			++n;
		}
	}
	// mean sliding window over full vector
	for ( ptrdiff_t i = 0; i < dst.len; ++i )
	{
		// update sum and n for src[i]
		if ( i > 0 )
		{
			ptrdiff_t prev = max2(i - hk - 1, 0);
			ptrdiff_t next = min2(i + hk, dst.len - 1);
			// neither are NA
			if ( !is_na(src[prev]) && !is_na(src[next]) )
			{
				sum = sum - src[prev] + src[next];
			}
			// next is NA
			else if ( !is_na(src[prev]) )
			{
				sum -= src[prev];
				--n;
			}
			// prev is NA
			else if ( !is_na(src[next]) )
			{
				sum += src[next];
				++n;
			}
		}
		// assign mean to dst[i]
		if ( n > 0 )
			dst[i] = sum / n;
		else
			dst[i] = na_value<T>();
	}
}

// Convolution filter
template<Num T, Vec V, Vec W>
void filt1_conv(vec<T> dst, V src, W w) noexcept
{
	dst.assign(convolve<T>(src, w));
}

#endif // CARDINAL_CORE_SIGNAL
