#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include "prelude.h"

//// Quicksort and Quickselect
//----------------------------

#define SWAP(x, y, T) do { T swap = x; x = y; y = swap; } while (false)
#define LINEAR_THRESHOLD 8

// select a pivot and partition x around the pivot such that
// * partitions the indices of x via out_index
// * all items left of pivot are <= pivot
// * all items right of pivot are >= pivot
// * incomparables sort last/highest (NA >> Inf)
// returns: pivot index
template<typename T, typename Index>
Index partition(
	const T * x,
	const Index lo, // index of first item to consider in x
	const Index hi, // index of last item to consider in x
	Index * out_index)
{
	// we get item k via x[at[k]]
	Index * at = out_index;
	// find pivot by median of 1st/mid/last
	Index pivot = static_cast<Index>((lo + hi) / 2);
	if ( LESSER(x[at[pivot]], x[at[lo]]) )
		SWAP(at[pivot], at[lo], Index);
	if ( GREATER(x[at[pivot]], x[at[hi]]) )
	{
		SWAP(at[pivot], at[hi], Index);
		if ( LESSER(x[at[pivot]], x[at[lo]]) )
			SWAP(at[pivot], at[lo], Index);
	}
	// lo and hi are now partitioned so skip them
	Index i = static_cast<Index>(lo + 1);
	Index j = static_cast<Index>(hi - 1);
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
// * sorts indices of x via out_index
// * sorts indices of elements within range
// * incomparables rank last/highest (NA >> Inf)
template<typename T, typename Index>
void quick_order(
	const vctr<T> x, 
	const slice range,
	Index * out_index,
	const bool init_out_index = false)
{
	// get the length of the slice
	size_t n = range.size();
	if ( n == 0 )
		return;
	// fill out_index with sequential indices
	if ( init_out_index )
		fill_buffer<Index>(out_index, n, 0, x.stride);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// initialize the stack
	size_t stack_size = 2 * std::ceil(std::log2(n) + 1);
	Index * stack = R_Calloc(stack_size, Index);
	Index top = -1;
	Index lo = static_cast<Index>(range.begin);
	Index hi = static_cast<Index>(range.end - 1);
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
		Index pivot = partition<T,Index>(x.ptr, lo, hi, at);
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
void quick_order(const vctr<T> x, Index * out_index)
{
	quick_order<T,Index>(x, x.all_elements(), out_index, true);
}

// find the k-th ranked item of an array x
// * partially sorts the indices of x via out_index
// * incomparables rank last/highest (NA >> Inf)
// returns: index of k-th item
template<typename T, typename Rank, typename Index>
T quick_select(
	const vctr<T> x,
	const slice range,
	const Rank k,
	Index * out_index,
	const bool init_out_index = false)
{
	// fill out_index with sequential indices
	if ( init_out_index )
		fill_buffer<Index>(out_index, range.size(), 0, x.stride);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// recursively partition the array
	Index lo = static_cast<Index>(range.begin);
	Index hi = static_cast<Index>(range.end - 1);
	do {
		if ( lo == hi )
			return x.ptr[at[lo]];
		Index pivot = partition<T,Index>(x.ptr, lo, hi, at);
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
// * incomparables rank last/highest (NA >> Inf)
// * returns the item values via out_values
template<typename T, typename Rank, typename Index>
void quick_select(const vctr<T> x, const vctr<Rank> k, T * out_values)
{
	// set up working index buffer
	Index * work_index = R_Calloc(x.len, Index);
	fill_buffer<Index>(work_index, x.len, 0, x.stride);
	// loop through k's	
	for ( size_t i = 0; i < k.len; ++i )
	{
		ptrdiff_t slen = static_cast<ptrdiff_t>(x.len);
		if ( i == 0 )
			out_values[0] = quick_select<T,Index,Rank>(
				x, slice{0, slen}, k.at(0), work_index);
		else if ( k.at(i) > k.at(i - 1) )
			out_values[i] = quick_select<T,Index,Rank>(
				x, slice{k.at(i - 1) + 1, slen}, k.at(i), work_index);
		else if ( k.at(i) < k.at(i - 1) )
			out_values[i] = quick_select<T,Index,Rank>(
				x, slice{0, k.at(i - 1)}, k.at(i), work_index);
		else
			out_values[i] = out_values[i - 1];
	}
	R_Free(work_index);
}

//// Median and MAD
//-----------------

// computes median of array x
// * incomparables are ignored/removed
// returns: the median
template<typename T, typename Index>
double quick_median(const vctr<T> x)
{
	// initialize result
	double median = mkIncomparable<double>();
	if ( x.len == 0 )
		return median;
	// set up working index buffer
	Index * work_index = R_Calloc(x.len, Index);
	fill_buffer<Index>(work_index, x.len, 0, x.stride);
	// find number of comparable items
	Index n = 0;
	for ( size_t i = 0; i < x.len; ++i )
	{
		if ( !isIncomparable(x.at(i)) )
			++n;
	}
	// compute median
	Index k = n / 2;
	if ( x.len % 2 == 0 )
	{
		ptrdiff_t len = static_cast<ptrdiff_t>(x.len);
		double m1 = quick_select<T,Index,Index>(
			x, slice{0, len}, k - 1, work_index);
		double m2 = quick_select<T,Index,Index>(
			x, slice{k, len}, k, work_index);
		median = 0.5 * (m1 + m2);
	}
	else
		median = quick_select<T,Index,Index>(
			x, x.all_elements(), k, work_index);
	R_Free(work_index);
	return median;
}

// computes MAD (Median Absolute Deviation) of array x
// * incomparables are ignored/removed
// * default scale is chosen so SD ~= MAD for if x ~ Normal
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
	for ( size_t i = 0; i < x.len; ++i )
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

// binary search for query x in data array
// * data must be sorted in non-decreasing order
// * differences <= tolerance are considered matches
// returns: index of match
template<typename T, typename Index>
Index binary_search(
	const T x,          // query
	const vctr<T> data, // data to search for query
	const double tolerance = DBL_EPSILON, 
	const bool relative = false, 
	const bool nearest = false,
	const Index nomatch = -1)
{
	if ( data.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = static_cast<Index>(data.len - 1);
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double d_mid = diff(data.at(mid), x, relative);
		if ( d_mid < 0 )
			lo = mid + 1;
		else if ( d_mid > 0 )
			hi = mid - 1;
		else
			return mid * data.stride;
	}
	double d_lo = std::fabs(diff(data.at(lo), x, relative));
	double d_hi = std::fabs(diff(data.at(hi), x, relative));
	if ( d_lo <= d_hi && (nearest || d_lo <= tolerance) )
		return lo;
	if ( d_hi <= d_lo && (nearest || d_hi <= tolerance) )
		return hi;
	return nomatch;
}

// binary search for multiple queries x in data array
// * data must be sorted in non-decreasing order
// * differences <= tolerance are considered matches
// * returns matches via out_index
template<typename T, typename Index>
void binary_search(
	const vctr<T> x, 
	const vctr<T> data,
	Index * out_index,
	const double tolerance = DBL_EPSILON, 
	const bool relative = false,
	const bool nearest = false, 
	const Index nomatch = -1)
{
	for ( size_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			out_index[i] = nomatch;
		else
		{
			out_index[i] = binary_search<T,Index>(
				x.at(i), 
				data, 
				tolerance, 
				relative, 
				nearest, 
				nomatch);
		}
	}
}

#endif // CARDINAL_CORE_SEARCH
