#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include <bit>
#include "core.h"
#include "order.h"

//// Binary search
//-----------------

// Binary search for query in x
// - Values of x MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<typename Index, typename T>
Index binary_search(
	const T query,
	const vec<T> x,
	const double tolerance = 0,
	const bool relative = false,
	const bool nearest = false,
	const Index nomatch = na_value<Index>())
{
	if ( x.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = x.len - 1;
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double dx = diff(x[mid], query, relative);
		if ( dx < 0 )
			lo = mid + 1;
		else if ( dx > 0 )
			hi = mid - 1;
		else
			return mid;
	}
	double dlo = std::fabs(diff(x[lo], query, relative));
	double dhi = std::fabs(diff(x[hi], query, relative));
	if ( dlo <= dhi && (nearest || dlo <= tolerance) )
		return lo;
	if ( dhi <= dlo && (nearest || dhi <= tolerance) )
		return hi;
	return nomatch;
}

// Binary search for multiple queries in x
// - Values of x MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<typename Index, typename T>
void binary_search(
	vec<Index> index,
	const vec<T> query,
	const vec<T> x,
	const double tolerance = 0,
	const bool relative = false,
	const bool nearest = false,
	const Index nomatch = na_value<Index>())
{
	for ( ptrdiff_t i = 0; i < query.len; ++i )
	{
		if ( is_na(query[i]) )
			index[i] = nomatch;
		else
		{
			index[i] = binary_search(
				query[i], 
				x, 
				tolerance, 
				relative, 
				nearest, 
				nomatch);
		}
	}
}

#endif // CARDINAL_CORE_SEARCH
