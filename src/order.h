#ifndef CARDINAL_CORE_ORDER
#define CARDINAL_CORE_ORDER

#include <bit>
#include <memory>
#include "core.h"

#define SMALL_SORT_THRESHOLD 8

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
	if ( is_na(x) && is_na(ref) )
		return 0.0;       // NAs sort equivalently
	else if ( is_na(x) )
		return +INFINITY; // NAs sort last so (x - ref) => +Inf
	else if ( is_na(ref) )
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
// Sorting and selection routines

// Swap items at indices i and j
// - Returns the vec for convenience
template<typename Index, typename T>
void swap(
	vec<T> x,
	const Index i,
	const Index j)
{
	T tmp = x[i];
	x[i] = x[j];
	x[j] = tmp;
}

// Select a pivot and partition x around the pivot such that
// - Caller MUST initialize out_index with valid indices of x
// - All items left of pivot are <= pivot
// - All items right of pivot are >= pivot
// - Incomparables sort last/highest (NA >> Inf)
// returns: pivot index
template<typename Index, typename T>
Index partition(
	vec<Index> index,
	const vec<T> x,
	const Index lo, // index of first item to consider in x
	const Index hi) // index of last item to consider in x
{
	// checks (debug only)
	assert(lo <= hi);
	assert(0 <= lo && lo < index.len);
	assert(0 <= hi && hi < index.len);
	// find pivot by median of 1st/mid/last
	Index pivot = (lo + hi) / 2;
	if ( LESSER(x[index[pivot]], x[index[lo]]) )
		swap(index, pivot, lo);
	if ( GREATER(x[index[pivot]], x[index[hi]]) ) {
		swap(index, pivot, hi);
		if ( LESSER(x[index[pivot]], x[index[lo]]) )
			swap(index, pivot, lo);
	}
	// lo and hi are now partitioned so skip them
	Index i = lo + 1;
	Index j = hi - 1;
	// use Hoare's partition method 
	do {
		// find next item less than pivot
		while ( LESSER(x[index[i]], x[index[pivot]]) ) ++i;
		// find next item greater than pivot
		while ( GREATER(x[index[j]], x[index[pivot]]) ) --j;
		// swap items (only if pointers haven't crossed)
		if ( i < j && NOT_EQUAL(x[index[i]], x[index[j]]) )
		{
			swap(index, i, j);
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
	vec<Index> index,
	const vec<T> x, 
	const bounds b)
{
	// checks (debug only)
	assert(0 <= b.start && b.start < index.len);
	assert(0 <= b.stop && b.stop <= index.len);
	assert(b.width() >= 0);
	// initialize the stack
	int stack_n = 2; // lo, hi
	int stack_size = stack_n * std::bit_width(static_cast<size_t>(b.width()));
	auto stack = std::make_unique<Index[]>(stack_size);
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
		if ( hi - lo < SMALL_SORT_THRESHOLD )
		{
			// use insertion sort for small subarrays
			for ( Index i = lo + 1; i <= hi; ++i )
			{
				Index j = i;
				while ( j > lo && LESSER(x[index[j]], x[index[j - 1]]) )
				{
					swap(index, j, j - 1);
					--j;
				}
			}
			// skip to next subarray
			continue;
		}
		Index pivot = partition(index, x, lo, hi);
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

template<typename Index, typename T>
void quick_order(
	vec<Index> index,
	const vec<T> x)
{
	quick_order(index, x, index.all_elements());
}

// Find the k-th ranked item of an array x
// - Caller SHOULD initialize out_index with valid indices of x
// - Partially sorts indices of x such that x[out_index[k]] is a pivot
// - Incomparables rank last/highest (NA >> Inf)
// returns: value of k-th item
template<typename Index, typename Rank, typename T>
T quick_select(
	vec<Index> index,
	const vec<T> x,
	const Rank k,
	const bounds b)
{
	// checks (debug only)
	assert(0 <= k && k < index.len);
	assert(0 <= b.start && b.start < index.len);
	assert(0 <= b.stop && b.stop <= index.len);
	assert(b.width() >= 0);
	// recursively partition the array
	Index lo = b.start;
	Index hi = b.stop - 1;
	do {
		if ( lo == hi )
			return x[index[lo]];
		Index pivot = partition(index, x, lo, hi);
		// return k-th element or partition again
		if ( k == pivot )
			return x[index[k]];
		else if ( k < pivot )
			hi = pivot - 1;
		else
			lo = pivot + 1;
	}
	while (true);
}

template<typename Index, typename Rank, typename T>
T quick_select(
	vec<Index> index,
	const vec<T> x,
	const Rank k)
{
	return quick_select(index, x, k, index.all_elements());
}

//// Median and MAD
//-----------------

// Computes median of array x
// - Incomparables are ignored/removed
// returns: the median
template<typename Index = ptrdiff_t, typename T>
double quick_median(const vec<T> x)
{
	// initialize result
	double median = na_value<double>();
	if ( x.len == 0 )
		return median;
	// set up working index buffer
	local_vec<Index> index{x.len};
	index.seqfill();
	// find count of comparable items
	Index n = 0;
	for ( Index i = 0; i < x.len; ++i )
	{
		if ( !is_na(x[i]) )
			++n;
	}
	// compute median
	Index k = n / 2;
	if ( x.len % 2 == 0 )
	{
		double m1 = quick_select(index.borrow(), x, k - 1);
		double m2 = quick_select(index.borrow(), x, k);
		return 0.5 * (m1 + m2);
	}
	else
	{
		return quick_select(index.borrow(), x, k);
	}
}

// Computes MAD (Median Absolute Deviation) of array x
// - Incomparables are ignored/removed
// - Default scale is chosen so SD ~= MAD for if x ~ Normal
// returns: the MAD
template<typename Index = ptrdiff_t, typename T>
double quick_mad(const vec<T> x, double center, double constant = 1.4826)
{
	if ( x.len == 0 )
		return na_value<double>();
	// compute absolute deviations
	local_vec<double> devs{x.len};
	for ( Index i = 0; i < x.len; ++i )
	{
		if ( is_na(x[i]) )
			devs[i] = na_value<double>();
		else
			devs[i] = std::fabs(x[i] - center);
	}
	return constant * quick_median<Index>(devs);
}

#endif // CARDINAL_CORE_ORDER
