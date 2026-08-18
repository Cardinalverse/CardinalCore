#ifndef CARDINAL_CORE_ORDER
#define CARDINAL_CORE_ORDER

#include <bit>
#include <memory>
#include "core.h"

#define SMALL_SORT_THRESHOLD 16

//// Order
//---------
// Traits supporting ordering

// An Ord Vec can be re-ordered and sorted
// - MUST have compare(i, j) -> 0 if o[i] == o[j]
// - MUST have compare(i, j) < 0 if o[i] < o[j]
// - MUST have compare(i, j) > 0 if o[i] < o[j]
// - SHOULD treat incomparables as sorting last and equivalently
template<class O>
concept Ord = Vec<O> &&
	requires(std::remove_cvref_t<O>& o, ptrdiff_t i, ptrdiff_t j)
	{
		{ o.compare(i, j) } -> Num;
		{ o.swap(i, j) };
	};

// An Ord Vec that re-orders indices instead of values
template<Vec V, Ord Index, Num T = double>
struct vec_ordered : vec_indexed<V,Index,T>
{
	T compare(const ptrdiff_t i, const ptrdiff_t j) const noexcept
	{
		T lhs = (*this)[i];
		T rhs = (*this)[j];
		if ( !is_na(lhs) && !is_na(rhs) )
			return lhs - rhs ? lhs - rhs : this->index[i] - this->index[j];
		else
			return is_na(lhs) - is_na(rhs);
	}

	void swap(const ptrdiff_t i, const ptrdiff_t j) noexcept
	{
		return this->index.swap(i, j);
	}
};

// A safe recusion depth limit for a balanced binary tree
template<Num T>
size_t max_depth(T n) {
	return 1 + std::bit_width(static_cast<size_t>(n));
}

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
	const bool relative = false) noexcept
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

//// Quicksort and Quickselect
//----------------------------
// Sorting and selection routines

// Select a pivot and partition range [x[lo]], x[hi]]
// - All x[i] left of pivot are <= pivot
// - All x[i] right of pivot are >= pivot
// - Missing/incomparable values sort last/highest (NA >> Inf)
// - Returns the pivot
template<Ord V>
ptrdiff_t partition(V x, const ptrdiff_t lo, const ptrdiff_t hi) noexcept
{
	// invariants
	assert(lo <= hi);
	assert(0 <= lo && lo < x.ssize());
	assert(0 <= hi && hi < x.ssize());
	// find pivot by median of 1st/mid/last
	ptrdiff_t pivot = (lo + hi) / 2;
	if ( x.compare(pivot, lo) < 0 )
		x.swap(pivot, lo);
	if ( x.compare(pivot, hi) > 0 ) {
		x.swap(pivot, hi);
		if ( x.compare(pivot, lo) < 0 )
			x.swap(pivot, lo);
	}
	// lo and hi are now partitioned so skip them
	ptrdiff_t i = lo + 1;
	ptrdiff_t j = hi - 1;
	// use Hoare's partition method
	do {
		// find next item not less than pivot
		while ( x.compare(i, pivot) < 0 ) ++i;
		// find next item not greater than pivot
		while ( x.compare(j, pivot) > 0 ) --j;
		// swap items (only if pointers haven't crossed)
		if ( i < j )
		{
			// swap pivot (if necessary)
			if ( pivot == i )
				pivot = j;
			else if ( pivot == j )
				pivot = i;
			// swap index at i and j
			x.swap(i, j);
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

// Sort an vector x
// - Sorts half-open interval [x[b.start]], x[b.stop])
// - Missing/incomparable values sort last/highest (NA >> Inf)
template<Ord V>
void qsort(V x, const bounds b)
{
	// invariants
	assert(b.start <= b.stop);
	assert(0 <= b.start && b.start < x.ssize());
	assert(0 <= b.stop && b.stop <= x.ssize());
	// initialize stack
	ptrdiff_t top = -1;
	struct frame { ptrdiff_t lo, hi; };
	auto stack = std::make_unique<frame[]>(max_depth(b.stop - b.start));
	stack[++top] = {b.start, b.stop - 1};
	// recursively partition the array
	while ( top >= 0 )
	{
		// pop and partition current subarray
		frame cur = stack[top--];
		if ( cur.hi - cur.lo < SMALL_SORT_THRESHOLD )
		{
			// use insertion sort for small subarrays
			for ( ptrdiff_t i = cur.lo + 1; i <= cur.hi; ++i )
			{
				ptrdiff_t j = i;
				while ( j > cur.lo && x.compare(j, j - 1) < 0 )
				{
					x.swap(j, j - 1);
					--j;
				}
			}
			// skip to next subarray
			continue;
		}
		ptrdiff_t pivot = partition(x, cur.lo, cur.hi);
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

template<Ord V>
void qsort(V x) { 
	qsort(x, {0, x.ssize()});
}

template<Ord Index, Vec V>
void qsort_index(Index index, const V x, const bounds b) {
	qsort(vec_ordered<V,Index>{x, index}, b);
}

template<Ord Index, Vec V>
void qsort_index(Index index, const V x) {
	qsort_index(index, x, {0, index.ssize()});
}

// Select k-th ranked element in vector x
// - Quickselect on half-open interval [x[b.start], x[b.stop])
// - Missing/incomparable values sort last/highest (NA >> Inf)
template<Ord V, Num Rank>
auto qselect(V x, const Rank k, const bounds b) noexcept
{
	// invariants
	assert(b.start <= b.stop);
	assert(0 <= b.start && b.start < x.ssize());
	assert(0 <= b.stop && b.stop <= x.ssize());
	assert(0 <= k && k < x.ssize());
	// recursively partition the array
	ptrdiff_t lo = b.start;
	ptrdiff_t hi = b.stop - 1;
	do {
		if ( lo == hi )
			return x[lo];
		ptrdiff_t pivot = partition(x, lo, hi);
		// return k-th element or partition again
		if ( k == pivot )
			return x[k];
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
	const V data,
	const Rank k) noexcept
{
	auto v = vec_ordered<V,vec<Index>>{data, index};
	return qselect(v, k, index.all_elements());
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
