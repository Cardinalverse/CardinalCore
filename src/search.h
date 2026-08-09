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
template<Num Index, Num T>
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
template<Num Index, Num T>
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

//// K-D search
//--------------

// KDTree for nearest neighbor and range searches
// - Owner is responsible for managing memory
// - Rows are observations and cols are features
// - Tree uses vectors giving indices of children
// - NA or negative indices indicate leaf nodes
template<Num Index, Num T>
struct kdtree
{
	mat<T> data;
	vec<Index> left;
	vec<Index> right;
	Index root = na_value<Index>();
	bool built = false;

	ptrdiff_t ssize() const noexcept
	{
		return data.nrows;
	}

	ptrdiff_t build() noexcept
	{
		// invariants
		assert(data.nrows == left.len);
		assert(data.nrows == right.len);
		if ( data.ssize() <= 0 )
			return root;
		// initialize indices
		local_vec<Index> index{data.nrows};
		index.seqfill(0);
		left.fill(na_value<Index>());
		right.fill(na_value<Index>());
		// find root from median of first dim
		vec<T> xcur = data.col(0);
		qsort_index(index, xcur);
		ptrdiff_t mid = data.nrows / 2;
		// handle duplicates and update root
		while ( mid > 0 && xcur[index[mid - 1]] == xcur[index[mid]] )
			--mid;
		root = index[mid];
		// initialize stack
		ptrdiff_t top = -1;
		struct frame { ptrdiff_t parent, depth, start, stop; };
		size_t max_depth = std::bit_width(2 * static_cast<size_t>(data.nrows));
		auto stack = std::make_unique<frame[]>(max_depth);
		// push initial left span to stack
		if ( mid > 0 ) {
			stack[++top] = {
				.parent = root,
				.depth = 1,
				.start = 0,
				.stop = mid,
			};
		}
		// push initial right span to stack
		if ( mid + 1 < data.nrows ) {
			stack[++top] = {
				.parent = root,
				.depth = 1,
				.start = mid + 1,
				.stop = data.nrows,
			};
		}
		// recursively build the tree
		while ( top >= 0 )
		{
			// pop stack
			frame cur = stack[top--];
			xcur = data.col(cur.depth % data.ncols);
			// find median of current dim within span
			if ( data.ncols > 1 )
				qsort_index(index, xcur, {cur.start, cur.stop});
			mid = (cur.start + cur.stop) / 2;
			while ( mid > 0 && xcur[index[mid - 1]] == xcur[index[mid]] )
				--mid;
			// insert child under parent
			ptrdiff_t pcol = (cur.depth - 1) % data.ncols;
			if ( data[{index[mid], pcol}] < data[{cur.parent, pcol}] )
				left[cur.parent] = index[mid];
			else
				right[cur.parent] = index[mid];
			// push left span
			if ( mid > cur.start ) {
				stack[++top] = {
					.parent = index[mid],
					.depth = cur.depth + 1,
					.start = cur.start,
					.stop = mid,
				};
			}
			// push right span
			if ( mid + 1 < cur.stop ) {
				stack[++top] = {
					.parent = index[mid],
					.depth = cur.depth + 1,
					.start = mid + 1,
					.stop = cur.stop,
				};
			}
		}
		built = true;
		return root;
	}

	#ifdef USING_R
	static kdtree<Index,T> from(
		SEXP data,
		SEXP left,
		SEXP right)
	{
		return {
			.data = mat<T>::from(data),
			.left = vec<Index>::from(left),
			.right = vec<Index>::from(right),
		};
	}
	#endif // USING_R
};

#endif // CARDINAL_CORE_SEARCH
