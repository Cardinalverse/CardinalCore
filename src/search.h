#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include <memory>
#include "core.h"
#include "order.h"

//// Utility
//-----------
// Search utilities

// Does query neighbor ref within some tolerance(s)?
template<Num T, Vec Tol, Vec Rel>
bool near(
	const vec<T> query,
	const vec<T> ref,
	const Tol tolerance,
	const Rel relative) noexcept
{
	assert(query.ssize() == ref.ssize());
	assert(query.ssize() == tolerance.ssize());
	assert(query.ssize() == relative.ssize());
	for ( ptrdiff_t i = 0; i < query.ssize(); ++i )
	{
		if ( std::fabs(diff(query[i], ref[i], relative[i])) > tolerance[i] )
			return false;
	}
	return true;
}

//// Binary search
//-----------------

// Binary search for query in red
// - Values of ref MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<Num Index, Num T>
Index binary_search(
	const T query,
	const vec<T> ref,
	const double tolerance = 0,
	const bool relative = false,
	const bool nearest = false,
	const Index nomatch = na_value<Index>())
{
	if ( ref.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = ref.len - 1;
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double dx = diff(query, ref[mid], relative);
		if ( dx > 0 )
			lo = mid + 1;
		else if ( dx < 0 )
			hi = mid - 1;
		else
			return mid;
	}
	double dlo = std::fabs(diff(query, ref[lo], relative));
	double dhi = std::fabs(diff(query, ref[hi], relative));
	if ( dlo <= dhi && (nearest || dlo <= tolerance) )
		return lo;
	if ( dhi <= dlo && (nearest || dhi <= tolerance) )
		return hi;
	return nomatch;
}

// Binary search for multiple queries in ref
// - Values of ref MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<Num Index, Num T>
void binary_search(
	vec<Index> index,
	const vec<T> query,
	const vec<T> ref,
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
				ref, 
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

	ptrdiff_t ssize() const noexcept
	{
		return data.nrows;
	}

	bool has_left(ptrdiff_t node) const noexcept
	{
		return 0 <= left[node] && left[node] < data.nrows;
	}

	bool has_right(ptrdiff_t node) const noexcept
	{
		return 0 <= right[node] && right[node] < data.nrows;
	}

	// Build the tree and return the index of the root node
	ptrdiff_t build()
	{
		// invariants
		assert(left.len == data.nrows);
		assert(right.len == data.nrows);
		if ( data.ssize() <= 0 )
			return root;
		// initialize indices
		local_vec<Index> index{data.nrows};
		index.seqfill(0);
		left.fill(na_value<Index>());
		right.fill(na_value<Index>());
		// find root from median of first dim
		vec<T> ref = data.col(0);
		qsort_index(index.borrow(), ref);
		ptrdiff_t mid = data.nrows / 2;
		// handle duplicates and update root
		while ( mid > 0 && ref[index[mid - 1]] == ref[index[mid]] )
			--mid;
		root = index[mid];
		// initialize stack
		ptrdiff_t top = -1;
		struct frame { ptrdiff_t parent, depth, start, stop; };
		auto stack = std::make_unique<frame[]>(max_depth(data.nrows));
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
			// get current dim
			ref = data.col(cur.depth % data.ncols);
			// find median of current dim within current span
			if ( data.ncols > 1 )
				qsort_index(index.borrow(), ref, {cur.start, cur.stop});
			mid = (cur.start + cur.stop) / 2;
			while ( mid > 0 && ref[index[mid - 1]] == ref[index[mid]] )
				--mid;
			// insert child under parent
			ptrdiff_t iprev = (cur.depth - 1) % data.ncols;
			if ( data[{index[mid], iprev}] < data[{cur.parent, iprev}] )
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
		return root;
	}

	// // Search for nearest neighbors of query
	// // - Number of neighbors (k) given by hits.len
	// // - Fill hits with indices up to hits.len
	// // - Fill dists with distances up to hits.len
	// template<Num Dist>
	// size_t knn_search(
	// 	vec<Index> hits,
	// 	vec<Dist> dists,
	// 	const vec<T> query) const
	// {
	// 	// invariants
	// 	assert(query.len == data.ncols);
	// 	assert(hits.len == dists.len)
	// 	// initialize stack
	// 	ptrdiff_t top = -1;
	// 	struct frame { ptrdiff_t node, depth; };
	// 	auto stack = std::make_unique<frame[]>(max_depth(data.nrows));
	// 	stack[++top] = { root, 0 };
	// 	// initialize hits
	// 	size_t count = 0;
	// 	hits.fill(na_value<Index>());
	// 	// recursively search tree
	// 	while ( top >= 0 )
	// 	{
	// 		// pop node
	// 		frame cur = stack[top--];
	// 		ptrdiff_t i = cur.depth % data.ncols;
	// 		double ds = diff(query[i], data[{cur.node, i}], relative[i]);
	// 		double du = std::fabs(ds);
	// 		// search left subtree?
	// 		if ( has_left(cur.node) && (ds < 0 || du <= tolerance[i]) )
	// 			stack[++top] = { left[cur.node], cur.depth + 1 };
	// 		// search right subtree?
	// 		if ( has_right(cur.node) && (ds > 0 || du <= tolerance[i]) )
	// 			stack[++top] = { right[cur.node], cur.depth + 1 };
	// 		// is this a hit?
	// 		if ( near(query, data.row(cur.node), tolerance, relative) )
	// 		{
	// 			if ( count < hits.len ) {
	// 				hits[count] = cur.node;
	// 				// sort into index order
	// 				ptrdiff_t j = count;
	// 				while ( j > 0 )
	// 				{
	// 					if ( hits[j] < hits[j - 1] )
	// 						hits.swap(j, j - 1);
	// 					--j;
	// 				}
	// 			}
	// 			++count;
	// 		}
	// 	}
	// 	return count;
	// }

	// Search for points within tolerance(s) of query
	// - Both tolerance and relative are per-dimension
	// - Fill hits with indices up to hits.len
	// - Return the count of hits
	template<Vec Tol, Vec Rel>
	size_t range_search(
		vec<Index> hits,
		const vec<T> query,
		const Tol tolerance,
		const Rel relative) const
	{
		// invariants
		assert(query.len == data.ncols);
		assert(tolerance.ssize() == data.ncols);
		assert(relative.ssize() == data.ncols);
		// initialize stack
		ptrdiff_t top = -1;
		struct frame { ptrdiff_t node, depth; };
		auto stack = std::make_unique<frame[]>(max_depth(data.nrows));
		stack[++top] = { root, 0 };
		// initialize hits
		size_t count = 0;
		hits.fill(na_value<Index>());
		// recursively search tree
		while ( top >= 0 )
		{
			// pop node
			frame cur = stack[top--];
			ptrdiff_t i = cur.depth % data.ncols;
			double ds = diff(query[i], data[{cur.node, i}], relative[i]);
			double du = std::fabs(ds);
			// search left subtree?
			if ( has_left(cur.node) && (ds < 0 || du <= tolerance[i]) )
				stack[++top] = { left[cur.node], cur.depth + 1 };
			// search right subtree?
			if ( has_right(cur.node) && (ds > 0 || du <= tolerance[i]) )
				stack[++top] = { right[cur.node], cur.depth + 1 };
			// is this a hit?
			if ( near(query, data.row(cur.node), tolerance, relative) )
			{
				if ( count < hits.len ) {
					hits[count] = cur.node;
					// sort into index order
					ptrdiff_t j = count;
					while ( j > 0 )
					{
						if ( hits[j] < hits[j - 1] )
							hits.swap(j, j - 1);
						--j;
					}
				}
				++count;
			}
		}
		return count;
	}

	template<Vec Tol, Vec Rel>
	size_t range_count(
		const vec<T> query,
		const Tol tolerance,
		const Rel relative) const
	{
		return range_search(
			vec<Index>{nullptr, 0, 0},
			query,
			tolerance,
			relative);
	}

	#ifdef USING_R
	static kdtree<Index,T> from(SEXP obj)
	{
		return {
			.data = mat<T>::from(VECTOR_ELT(obj, 0)),
			.left = vec<Index>::from(VECTOR_ELT(obj, 1)),
			.right = vec<Index>::from(VECTOR_ELT(obj, 2)),
			.root = *data_ptr<Index>(VECTOR_ELT(obj, 3)),
		};
	}
	#endif // USING_R
};

template<Num Index, Num T, Vec Tol, Vec Rel>
struct range_counts
{
	kdtree<Index,T> searcher;
	vec<Index> counts;
	mat<T> query;
	Tol tolerance;
	Rel relative;

	ptrdiff_t ssize() const
	{
		return query.nrows;
	}

	void operator()(bounds b)
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
		{
			counts[i] = searcher.range_count(
				query.row(i),
				tolerance,
				relative);
		}
	}
};

template<Num Index, Num T, Vec Tol, Vec Rel>
struct range_searches
{
	kdtree<Index,T> searcher;
	rag<Index,Index> hits;
	mat<T> query;
	Tol tolerance;
	Rel relative;

	ptrdiff_t ssize() const
	{
		return query.nrows;
	}

	void operator()(bounds b)
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
		{
			searcher.range_search(
				hits[i],
				query.row(i),
				tolerance,
				relative);
		}
	}
};

#endif // CARDINAL_CORE_SEARCH
