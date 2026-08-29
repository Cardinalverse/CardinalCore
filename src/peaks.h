#ifndef CARDINAL_CORE_PEAKS
#define CARDINAL_CORE_PEAKS

#include "core.h"

//// Peak detection
//------------------
// Find peaks in a signal

// Apply Callable to peaks in a signal
// - A peak is a local maximum among k points
// - Callable will be called as f(Index)
// - Returns the count of peaks
template<UnaryOp Callable, Vec V, Num Index>
Index peaks_apply(
	Callable f,
	const V x,
	const Index k = 5) noexcept
{
	if ( x.ssize() < k )
		return 0;
	Index halfWindow = k / 2;
	Index start = halfWindow;
	Index stop = x.ssize() - halfWindow;
	Index count = 0;
	for ( Index i = start; i < stop; ++i )
	{
		bool is_peak = true;
		// Peak must be > all points to left
		for ( Index j = i - halfWindow; j < i; ++j )
		{
			if ( x[j] >= x[i] )
			{
				is_peak = false;
				break;
			}
		}
		// Peak must be >= all points to right
		for ( Index j = i + 1; j <= i + halfWindow; ++j )
		{
			if ( x[j] > x[i] )
			{
				is_peak = false;
				break;
			}
		}
		// Process peak index
		if ( is_peak )
		{
			f(i);
			++count;
		}
	}
	return count;
}

// Count the number of peaks in a signal
template<Vec V, Num Index>
Index peaks_count(const V x, const Index k = 5) noexcept
{
	return peaks_apply(sink<Index>{}, x, k);
}

// Find the indices of peaks in a signal
// - Fills index with hits (up to index.len)
// - Returns count of hits (may be > index.len)
template<Vec V, Num Index>
Index peaks_find(
	vec<Index> index,
	const V x,
	const Index k = 5) noexcept
{
	return peaks_apply(sink<Index>{index}, x, k);
}

#endif // CARDINAL_CORE_PEAKS
