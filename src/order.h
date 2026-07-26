#ifndef CARDINAL_CORE_ORDER
#define CARDINAL_CORE_ORDER

#include <bit>
#include "core.h"

#define SWAP(x, y, T) do { T swap = x; x = y; y = swap; } while (false)
#define LINEAR_THRESHOLD 8

//// Comparison
//--------------
// Comparisons handling incomparables (NAs and NaNs)

// compute signed absolute or relative difference
// - safe to use with incomparables (NAs and NaNs)
// - incomparables sort last/highest (NA >> Inf)
// - incomparables sort equal to each other (NA == NA)
// returns: the difference
template<typename T>
double diff(
	const T x, 
	const T ref, 
	const bool relative = false)
{
	if ( isIncomparable(x) && isIncomparable(ref) )
		return 0.0;       // NAs sort equivalently
	else if ( isIncomparable(x) )
		return +INFINITY; // NAs sort last so (x - ref) => +Inf
	else if ( isIncomparable(ref) )
		return -INFINITY; // NAs sort last so (x - ref) => -Inf
	else
	{
		if ( relative )
			return static_cast<double>(x - ref) / ref;
		else
			return static_cast<double>(x - ref);
	}
}

#define LESSER(x, y) (diff((x), (y)) < 0)
#define GREATER(x, y) (diff((x), (y)) > 0)
#define LESSER_EQUAL(x, y) (diff((x), (y)) <= 0)
#define GREATER_EQUAL(x, y) (diff((x), (y)) >= 0)
#define EQUAL(x, y) (diff((x), (y)) == 0)
#define NOT_EQUAL(x, y) (diff((x), (y)) != 0)

//// Quicksort and Quickselect
//----------------------------

// Select a pivot and partition x around the pivot such that
// - Caller MUST initialize out_index with valid indices of x
// - All items left of pivot are <= pivot
// - All items right of pivot are >= pivot
// - Incomparables sort last/highest (NA >> Inf)
// returns: pivot index
template<typename Index, typename T>
Index partition(
	Index * out_index,
	const vec<T> x,
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

// Sort indices of an array x using Hoare's quicksort algorithm
// - Caller SHOULD initialize out_index with valid indices of x
// - Sorts indices of x such that x[out_index[i]] are sorted for i in slice
// - Incomparables rank last/highest (NA >> Inf)
template<typename Index, typename T>
void quick_order(
	Index * out_index,
	const vec<T> x, 
	const bounds b,
	const bool init_index = false)
{
	// check length of slice
	if ( b.len() <= 0 )
		return;
	// fill out_index with sequential indices
	if ( init_index )
		fill_buffer<Index>(out_index, x.len, 0, 1);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// initialize the stack
	int stack_n = 2; // lo, hi
	int stack_size = stack_n * std::bit_width(static_cast<size_t>(b.len()));
	Index * stack = SAFE_ALLOC(stack_size, Index);
	Index top = -1;
	Index lo = b.start;
	Index hi = b.stop - 1;
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
				while ( j > lo && LESSER(x[at[j]], x[at[j - 1]]) )
				{
					SWAP(at[j], at[j - 1], Index);
					--j;
				}
			}
			// skip to next subarray
			continue;
		}
		Index pivot = partition(at, x, lo, hi);
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
	SAFE_FREE(stack);
}

template<typename Index, typename T>
void quick_order(
	Index * out_index,
	const vec<T> x)
{
	quick_order(out_index, x, x.all_elements(), true);
}

// Find the k-th ranked item of an array x
// - Caller SHOULD initialize out_index with valid indices of x
// - Partially sorts indices of x such that x[out_index[k]] is a pivot
// - Incomparables rank last/highest (NA >> Inf)
// returns: value of k-th item
template<typename Index, typename Rank, typename T>
T quick_select(
	Index * out_index,
	const Rank k,
	const vec<T> x,
	const bounds b,
	const bool init_index = false)
{
	// fill out_index with sequential indices
	if ( init_index )
		fill_buffer<Index>(out_index, x.len, 0, 1);
	// we get item k via x[at[k]]
	Index * at = out_index;
	// recursively partition the array
	Index lo = b.start;
	Index hi = b.stop - 1;
	do {
		if ( lo == hi )
			return x[at[lo]];
		Index pivot = partition(at, x, lo, hi);
		// return k-th element or partition again
		if ( k == pivot )
			return x[at[k]];
		else if ( k < pivot )
			hi = pivot - 1;
		else
			lo = pivot + 1;
	}
	while (true);
}

// Find the k-th ranked items of an array for multiple k's
// - Caller SHOULD initialize out_index with valid indices of x
// - Partially sorts indices of x such that x[out_index[k]] is a pivot
// - Incomparables rank last/highest (NA >> Inf)
template<typename Rank, typename T>
void quick_select(
	T * out_values, 
	const vec<Rank> k,
	const vec<T> x)
{
	// set up working index buffer
	ptrdiff_t * index = SAFE_ALLOC(x.len, ptrdiff_t);
	fill_buffer<ptrdiff_t>(index, x.len, 0, 1);
	// loop through k's	
	for ( ptrdiff_t i = 0; i < k.len; ++i )
	{
		if ( i == 0 )
			out_values[0] = quick_select(
				index, k[0], x, {0, x.len});
		else if ( k[i] < k[i - 1] )
			out_values[i] = quick_select(
				index, k[i], x, {0, k[i - 1]});
		else if ( k[i] > k[i - 1] )
			out_values[i] = quick_select(
				index, k[i], x, {k[i - 1] + 1, x.len});
		else
			out_values[i] = out_values[i - 1];
	}
	SAFE_FREE(index);
}

//// Median and MAD
//-----------------

// Computes median of array x
// - Incomparables are ignored/removed
// returns: the median
template<typename T>
double quick_median(const vec<T> x)
{
	// initialize result
	double median = mkIncomparable<double>();
	if ( x.len == 0 )
		return median;
	// set up working index buffer
	ptrdiff_t * index = SAFE_ALLOC(x.len, ptrdiff_t);
	fill_buffer<ptrdiff_t>(index, x.len, 0, 1);
	// find count of comparable items
	ptrdiff_t count = 0;
	for ( ptrdiff_t i = 0; i < x.len; ++i )
	{
		if ( !isIncomparable(x[i]) )
			++count;
	}
	// compute median
	ptrdiff_t k = count / 2;
	if ( x.len % 2 == 0 )
	{
		double m1 = quick_select(index, k - 1, x, {0, x.len});
		double m2 = quick_select(index, k, x, {k, x.len});
		median = 0.5 * (m1 + m2);
	}
	else
		median = quick_select(index, k, x, x.all_elements());
	SAFE_FREE(index);
	return median;
}

// Computes MAD (Median Absolute Deviation) of array x
// - Incomparables are ignored/removed
// - Default scale is chosen so SD ~= MAD for if x ~ Normal
// returns: the MAD
template<typename T>
double quick_mad(const vec<T> x, double center, double scale = 1.4826)
{
	// initialize result
	double mad = mkIncomparable<double>();
	if ( x.len == 0 )
		return mad;
	// compute absolute deviations
	double * dev = SAFE_ALLOC(x.len, double);
	for ( ptrdiff_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x[i]) )
			dev[i] = mkIncomparable<double>();
		else
			dev[i] = std::fabs(x[i] - center);
	}
	mad = scale * quick_median(vec{dev, x.len, 1});
	SAFE_FREE(dev);
	return mad;
}

#endif // CARDINAL_CORE_ORDER
