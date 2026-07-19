#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include "prelude.h"

//// Quicksort and Quickselect
//----------------------------

#define SWAP(x, y, T) do { T swap = x; x = y; y = swap; } while (false)
#define LINEAR_THRESHOLD 8

// select a pivot and partition x around the pivot such that
// - partitions the indices of x via out_index
// - all items left of pivot are <= pivot
// - all items right of pivot are >= pivot
// - incomparables sort last/highest (NA >> Inf)
// returns: pivot index
template<typename T, typename Index>
Index partition(
	Index * out_index,
	const T * x,
	const Index lo, // index of first item to consider in x
	const Index hi) // index of last item to consider in x
{
	// we get item k via x[at[k]]
	Index * at = out_index;
	// find pivot by median of 1st/mid/last
	Index pivot = (lo + hi) / 2;
	if ( LESSER(x[at[pivot]], x[at[lo]]) )
		SWAP(at[pivot], at[lo], Index);
	if ( GREATER(x[at[pivot]], x[at[hi]]) )
	{
		SWAP(at[pivot], at[hi], Index);
		if ( LESSER(x[at[pivot]], x[at[lo]]) )
			SWAP(at[pivot], at[lo], Index);
	}
	// lo and hi are now partitioned so skip them
	Index i = lo + 1;
	Index j = hi - 1;
	// use Hoare's partition method 
	do {
		// find next item less than pivot
		while ( LESSER(x[at[i]], x[at[pivot]]) ) ++i;
		// find next item greater than pivot
		while ( GREATER(x[at[j]], x[at[pivot]]) ) --j;
		// swap items (only if pointers haven't crossed)
		if ( i < j && NOT_EQUAL(x[at[i]], x[at[j]]) )
		{
			SWAP(at[i], at[j], Index);
			if ( pivot == i )
				pivot = j;
			else if ( pivot == j )
				pivot = i;
		}
		// allow pointers to cross
		else if ( i == j )
		{
			++i;
			--j;
		}
		// account for ties
		else
		{
			if ( i != pivot )
				++i;
			if ( j != pivot )
				--j;
		}
	} while (i <= j);
	return pivot;
}

// sort an array x using Hoare's quicksort algorithm
// - sorts indices of x via out_index
// - sorts indices of elements within region
// - incomparables rank last/highest (NA >> Inf)
// - region.width must be nonnegative
template<typename T, typename Index>
void quick_order(
	Index * out_index,
	const vctr<T> x, 
	const slice region,
	const bool init_index = false)
{
	// get the length of the slice
	isize n = region.width();
	if ( n == 0 )
		return;
	// fill out_index with sequential indices
	if ( init_index )
		fill_buffer<Index>(out_index, n, 0, x.stride);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// initialize the stack
	isize stack_size = 2 * std::ceil(std::log2(n) + 1);
	Index * stack = R_Calloc(stack_size, Index);
	Index top = -1;
	Index lo = static_cast<Index>(region.start);
	Index hi = static_cast<Index>(region.stop - 1);
	stack[++top] = lo;
	stack[++top] = hi;
	// recursively partition the array
	while ( top >= 0 )
	{
		// pop and partition current subarray
		hi = stack[top--];
		lo = stack[top--];
		if ( hi - lo < LINEAR_THRESHOLD )
		{
			// use insertion sort for small subarrays
			for ( Index i = lo + 1; i <= hi; ++i )
			{
				Index j = i;
				while ( j > lo && LESSER(x.ptr[at[j]], x.ptr[at[j - 1]]) )
				{
					SWAP(at[j], at[j - 1], Index);
					--j;
				}
			}
			// skip to next subarray
			continue;
		}
		Index pivot = partition<T,Index>(at, x.ptr, lo, hi);
		// push larger subarray then smaller subarray
		if ( pivot - lo < hi - pivot )
		{
			// push higher subarray if non-empty
			if ( pivot + 1 < hi )
			{
				stack[++top] = pivot + 1;
				stack[++top] = hi;
			}
			// push lower subarray if non-empty
			if ( pivot - 1 > lo )
			{
				stack[++top] = lo;
				stack[++top] = pivot - 1;
			}
		}
		else {
			// push lower subarray if non-empty
			if ( pivot - 1 > lo )
			{
				stack[++top] = lo;
				stack[++top] = pivot - 1;
			}
			// push higher subarray if non-empty
			if ( pivot + 1 < hi )
			{
				stack[++top] = pivot + 1;
				stack[++top] = hi;
			}
		}
	}
}

template<typename T, typename Index>
void quick_order(Index * out_index, const vctr<T> x)
{
	quick_order<T,Index>(out_index, x, x.all_elements(), true);
}

// find the k-th ranked item of an array x
// - partially sorts the indices of x via out_index
// - incomparables rank last/highest (NA >> Inf)
// - region.width must be nonnegative
// returns: index of k-th item
template<typename T, typename Rank, typename Index>
T quick_select(
	Index * out_index,
	const vctr<T> x,
	const slice region,
	const Rank k,
	const bool init_index = false)
{
	// fill out_index with sequential indices
	if ( init_index )
		fill_buffer<Index>(out_index, region.width(), 0, x.stride);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// recursively partition the array
	Index lo = static_cast<Index>(region.start);
	Index hi = static_cast<Index>(region.stop - 1);
	do {
		if ( lo == hi )
			return x.ptr[at[lo]];
		Index pivot = partition<T,Index>(at, x.ptr, lo, hi);
		// return k-th element or partition again
		if ( k == pivot )
			return x.ptr[at[k]];
		else if ( k < pivot )
			hi = pivot - 1;
		else
			lo = pivot + 1;
	}
	while (true);
}

// find the k-th ranked items of an array for multiple k's
// - incomparables rank last/highest (NA >> Inf)
// - returns the item values via out_values
template<typename T, typename Rank, typename Index>
void quick_select(T * out_values, const vctr<T> x, const vctr<Rank> k)
{
	// set up working index buffer
	Index * index = R_Calloc(x.len, Index);
	fill_buffer<Index>(index, x.len, 0, x.stride);
	// loop through k's	
	for ( isize i = 0; i < k.len; ++i )
	{
		if ( i == 0 )
			out_values[0] = quick_select<T,Index,Rank>(
				 index, x, slice{0, x.len}, k.at(0));
		else if ( k.at(i) > k.at(i - 1) )
			out_values[i] = quick_select<T,Index,Rank>(
				index,  x, slice{k.at(i - 1) + 1, x.len}, k.at(i));
		else if ( k.at(i) < k.at(i - 1) )
			out_values[i] = quick_select<T,Index,Rank>(
				index, x, slice{0, k.at(i - 1)}, k.at(i));
		else
			out_values[i] = out_values[i - 1];
	}
	R_Free(index);
}

//// Median and MAD
//-----------------

// computes median of array x
// - incomparables are ignored/removed
// returns: the median
template<typename T, typename Index>
double quick_median(const vctr<T> x)
{
	// initialize result
	double median = mkIncomparable<double>();
	if ( x.len == 0 )
		return median;
	// set up working index buffer
	Index * index = R_Calloc(x.len, Index);
	fill_buffer<Index>(index, x.len, 0, x.stride);
	// find number of comparable items
	Index n = 0;
	for ( isize i = 0; i < x.len; ++i )
	{
		if ( !isIncomparable(x.at(i)) )
			++n;
	}
	// compute median
	Index k = n / 2;
	if ( x.len % 2 == 0 )
	{
		isize len = x.len;
		double m1 = quick_select<T,Index,Index>(
			index, x, slice{0, len}, k - 1);
		double m2 = quick_select<T,Index,Index>(
			index, x, slice{k, len}, k);
		median = 0.5 * (m1 + m2);
	}
	else
		median = quick_select<T,Index,Index>(
			index, x, x.all_elements(), k);
	R_Free(index);
	return median;
}

// computes MAD (Median Absolute Deviation) of array x
// - incomparables are ignored/removed
// - default scale is chosen so SD ~= MAD for if x ~ Normal
// returns: the MAD
template<typename T, typename Index>
double quick_mad(const vctr<T> x, double center, double scale = 1.4826)
{
	// initialize result
	double mad = mkIncomparable<double>();
	if ( x.len == 0 )
		return mad;
	// compute absolute deviations
	double * dev = R_Calloc(x.len, double);
	for ( isize i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			dev[i] = mkIncomparable<double>();
		else
			dev[i] = std::fabs(x.at(i) - center);
	}
	mad = scale * quick_median<double,Index>({dev, x.len, 1});
	R_Free(dev);
	return mad;
}

//// Binary search
//-----------------

// binary search for query in x
// - x must be sorted in non-decreasing order
// - differences <= tolerance are considered matches
// returns: index of match
template<typename T, typename Index>
Index binary_search(
	const T query, 
	const vctr<T> x, 
	const double tolerance = DBL_EPSILON, 
	const bool relative = false, 
	const bool nearest = false,
	const Index nomatch = -1)
{
	if ( x.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = static_cast<Index>(x.len - 1);
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double d_mid = diff(x.at(mid), query, relative);
		if ( d_mid < 0 )
			lo = mid + 1;
		else if ( d_mid > 0 )
			hi = mid - 1;
		else
			return mid * x.stride;
	}
	double d_lo = std::fabs(diff(x.at(lo), query, relative));
	double d_hi = std::fabs(diff(x.at(hi), query, relative));
	if ( d_lo <= d_hi && (nearest || d_lo <= tolerance) )
		return lo;
	if ( d_hi <= d_lo && (nearest || d_hi <= tolerance) )
		return hi;
	return nomatch;
}

// binary search for multiple queries x in data array
// - data must be sorted in non-decreasing order
// - differences <= tolerance are considered matches
// - returns matches via out_index
template<typename T, typename Index>
void binary_search(
	Index * out_index,
	const vctr<T> query, 
	const vctr<T> x,
	const double tolerance = DBL_EPSILON, 
	const bool relative = false,
	const bool nearest = false, 
	const Index nomatch = -1)
{
	for ( isize i = 0; i < query.len; ++i )
	{
		if ( isIncomparable(query.at(i)) )
			out_index[i] = nomatch;
		else
		{
			out_index[i] = binary_search<T,Index>(
				query.at(i), 
				x, 
				tolerance, 
				relative, 
				nearest, 
				nomatch);
		}
	}
}

#endif // CARDINAL_CORE_SEARCH
