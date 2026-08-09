#ifndef CARDINAL_CORE_ORDER
#define CARDINAL_CORE_ORDER

#include <bit>
#include <memory>
#include "core.h"

#define SMALL_SORT_THRESHOLD 16

//// Comparison
//--------------
// Comparisons handling NA/missing/incomparable

// Compute signed absolute or relative difference
// - Safe to use with incomparables (NAs and NaNs)
// - Incomparables sort last/highest (NA >> Inf)
// - Incomparables sort equal to each other (NA == NA)
template<Num T>
double diff(
	const T x, 
	const T ref, 
	const bool relative = false)
{
	if constexpr ( HasNA<T> )
	{
		if ( is_na(x) && is_na(ref) )
			return 0;         // NAs sort equivalently
		else if ( is_na(x) )
			return +INFINITY; // NAs sort last so (x - ref) => +Inf
		else if ( is_na(ref) )
			return -INFINITY; // NAs sort last so (x - ref) => -Inf
	}
	if ( relative )
		return static_cast<double>(x - ref) / ref;
	else
		return static_cast<double>(x - ref);
}

// Tests x[index[i]] < x[index[j]]
// - Ties broekn by index[i] < index[j]
// - NAs sort last (NA >> Inf)
template<Num Index, Vec V>
bool less_at(
	vec<Index> index,
	V x,
	ptrdiff_t i,
	ptrdiff_t j) noexcept
{
	Index ii = index[i];
	Index jj = index[j];
	double dx = diff(x[ii], x[jj]);
	return dx != 0 ? dx < 0 : ii < jj;
}

// Tests x[index[i]] > x[index[j]]
// - Ties broekn by index[i] > index[j]
// - NAs sort last (NA >> Inf)
template<Num Index, Vec V>
bool greater_at(
	vec<Index> index,
	V x,
	ptrdiff_t i,
	ptrdiff_t j) noexcept
{
	Index ii = index[i];
	Index jj = index[j];
	double dx = diff(x[ii], x[jj]);
	return dx != 0 ? dx > 0 : ii > jj;
}

//// Quicksort and Quickselect
//----------------------------
// Sorting and selection routines

// Select a pivot and partition range [x[index[lo]], x[index[hi]]]
// - All x[index[i]] left of pivot are <= pivot
// - All x[index[i]] right of pivot are >= pivot
// - Missing/incomparable values sort last/highest (NA >> Inf)
// - The index MUST contain valid indices of x
// - Returns the pivot
template<Num Index, Vec V>
ptrdiff_t partition_index(
	vec<Index> index,
	const V x,
	const ptrdiff_t lo, // index of first item to consider in x
	const ptrdiff_t hi) // index of last item to consider in x
{
	// invariants
	assert(lo <= hi);
	assert(0 <= lo && lo < index.len);
	assert(0 <= hi && hi < index.len);
	// find pivot by median of 1st/mid/last
	ptrdiff_t pivot = (lo + hi) / 2;
	if ( less_at(index, x, pivot, lo) )
		index.swap(pivot, lo);
	if ( greater_at(index, x, pivot, hi) ) {
		index.swap(pivot, hi);
		if ( less_at(index, x, pivot, lo) )
			index.swap(pivot, lo);
	}
	// lo and hi are now partitioned so skip them
	ptrdiff_t i = lo + 1;
	ptrdiff_t j = hi - 1;
	// use Hoare's partition method
	do {
		// find next item not less than pivot
		while ( less_at(index, x, i, pivot) ) ++i;
		// find next item not greater than pivot
		while ( greater_at(index, x, j, pivot) ) --j;
		// swap items (only if pointers haven't crossed)
		if ( i < j )
		{
			// swap pivot (if necessary)
			if ( pivot == i )
				pivot = j;
			else if ( pivot == j )
				pivot = i;
			// swap index at i and j
			index.swap(i, j);
		}
		else
		{
			// allow pointers to cross
			if ( i == j )
			{
				++i;
				--j;
			}
			// handle duplicate indices
			else
			{
				if ( i != pivot ) ++i;
				if ( j != pivot ) --j;
			}
		}
	} while (i <= j);
	return pivot;
}

// Sort indices of an array x
// - Sorts half-open interval [x[index[b.start]], x[index[b.stop]])
// - Missing/incomparable values sort last/highest (NA >> Inf)
// - The index MUST contain valid indices of x
template<Num Index, Vec V>
void qsort_index(
	vec<Index> index,
	const V x, 
	const bounds b)
{
	// invariants
	assert(b.start <= b.stop);
	assert(0 <= b.start && b.start < index.len);
	assert(0 <= b.stop && b.stop <= index.len);
	// initialize stack
	ptrdiff_t top = -1;
	struct frame { ptrdiff_t lo, hi; };
	size_t max_depth = std::bit_width(static_cast<size_t>(b.stop - b.start));
	auto stack = std::make_unique<frame[]>(max_depth);
	stack[++top] = {b.start, b.stop - 1};
	// recursively partition the array
	while ( top >= 0 )
	{
		// pop and partition current subarray
		frame cur = stack[top--];
		if ( cur.hi - cur.lo <= SMALL_SORT_THRESHOLD )
		{
			// use insertion sort for small subarrays
			for ( ptrdiff_t i = cur.lo + 1; i <= cur.hi; ++i )
			{
				ptrdiff_t j = i;
				while ( j > cur.lo && less_at(index, x, j, j - 1) )
				{
					index.swap(j, j - 1);
					--j;
				}
			}
			// skip to next subarray
			continue;
		}
		ptrdiff_t pivot = partition_index(index, x, cur.lo, cur.hi);
		// push larger subarray then smaller subarray
		if ( pivot - cur.lo < cur.hi - pivot )
		{
			// left < right => push right, then left
			if ( pivot + 1 < cur.hi )
				stack[++top] = {pivot + 1, cur.hi};
			if ( pivot - 1 > cur.lo )
				stack[++top] = {cur.lo, pivot - 1};
		}
		else
		{
			// left > right => push left, then right
			if ( pivot - 1 > cur.lo )
				stack[++top] = {cur.lo, pivot - 1};
			if ( pivot + 1 < cur.hi )
				stack[++top] = {pivot + 1, cur.hi};
		}
	}
}

template<Num Index, Vec V>
void qsort_index(vec<Index> index, const V x)
{
	qsort_index(index, x, index.all_elements());
}

// Select k-th ranked index in array x
// - Quickselect on half-open interval [x[index[b.start]], x[index[b.stop]])
// - Missing/incomparable values sort last/highest (NA >> Inf)
// - The index MUST contain valid indices of x
template<Num Index, Vec V, Num Rank>
auto qselect_index(
	vec<Index> index,
	const V x,
	const Rank k,
	const bounds b)
{
	// invariants
	assert(b.start <= b.stop);
	assert(0 <= b.start && b.start < index.len);
	assert(0 <= b.stop && b.stop <= index.len);
	assert(0 <= k && k < index.len);
	// recursively partition the array
	ptrdiff_t lo = b.start;
	ptrdiff_t hi = b.stop - 1;
	do {
		if ( lo == hi )
			return x[index[lo]];
		ptrdiff_t pivot = partition_index(index, x, lo, hi);
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

template<Num Index, Vec V, Num Rank>
auto qselect_index(
	vec<Index> index,
	const V x,
	const Rank k)
{
	return qselect_index(index, x, k, index.all_elements());
}

//// Median and MAD
//-----------------

// Computes median of array x
// - NA/missing/incomparable values are ignored
template<Vec V>
double qmedian(const V x)
{
	if ( x.ssize() == 0 )
		return na_value<double>();
	// initialize index
	local_vec<ptrdiff_t> index{x.ssize()};
	index.seqfill(0);
	// find count of non-NA items
	ptrdiff_t n = n_present(x);
	// compute median
	ptrdiff_t k = n / 2;
	if ( n % 2 == 0 )
	{
		double m1 = qselect_index(index.borrow(), x, k - 1);
		double m2 = qselect_index(index.borrow(), x, k);
		return 0.5 * (m1 + m2);
	}
	else
	{
		return qselect_index(index.borrow(), x, k);
	}
}

// Computes MAD (Median Absolute Deviation) of array x
// - Default constant is chosen so SD ~= MAD for if x ~ Normal
template<Vec V>
double qmad(const V x, double center, double constant = 1.4826)
{
	if ( x.ssize() == 0 )
		return na_value<double>();
	auto devs = ufunc<Abs>(x - rep{center, x.ssize()});
	return constant * qmedian(devs);
}

#endif // CARDINAL_CORE_ORDER
