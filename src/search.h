#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include <bit>
#include "core.h"
#include "order.h"

//// Binary search
//-----------------

// Binary search for query in x
// - Caller MUST guarantee values of x are non-decreasing
// - Differences <= tolerance are considered matches
// returns: index of match
template<typename Index, typename T>
Index binary_search(
	const T query, 
	const vec<T> x, 
	const double tolerance = DBL_EPSILON, 
	const bool relative = false, 
	const bool nearest = false,
	const Index nomatch = -1)
{
	if ( x.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = x.len - 1;
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double d_mid = diff(x[mid], query, relative);
		if ( d_mid < 0 )
			lo = mid + 1;
		else if ( d_mid > 0 )
			hi = mid - 1;
		else
			return mid;
	}
	double d_lo = std::fabs(diff(x[lo], query, relative));
	double d_hi = std::fabs(diff(x[hi], query, relative));
	if ( d_lo <= d_hi && (nearest || d_lo <= tolerance) )
		return lo;
	if ( d_hi <= d_lo && (nearest || d_hi <= tolerance) )
		return hi;
	return nomatch;
}

// Binary search for multiple queries in x
// - Caller MUST guarantee values of x are non-decreasing
// - Differences <= tolerance are considered matches
template<typename Index, typename T>
void binary_search(
	Index * out_index,
	const vec<T> query, 
	const vec<T> x,
	const double tolerance = DBL_EPSILON, 
	const bool relative = false,
	const bool nearest = false, 
	const Index nomatch = -1)
{
	for ( ptrdiff_t i = 0; i < query.len; ++i )
	{
		if ( isIncomparable(query[i]) )
			out_index[i] = nomatch;
		else
		{
			out_index[i] = binary_search(
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
