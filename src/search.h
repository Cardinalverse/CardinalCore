#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

#include <memory>
#include "core.h"
#include "kernels.h"
#include "order.h"
#include "dist.h"

//// Utility
//-----------
// Search utilities

enum Diff {
	Absolute, // Absolute diff
	RefLhs,   // Relative diff using Lhs as ref
	RefRhs,   // Relative diff using Rhs as ref
};

enum Ref {
	Query, // Use search 'query' as ref (for relative diff)
	Table, // Use search 'table' as ref (for relative diff)
};

// Compute signed absolute or relative difference
template<Diff Method = Absolute, Num L, Num R>
double diff(const L lhs, const R rhs) noexcept
{
	if ( is_na(lhs) || is_na(rhs) )
		return huge_positive_value<double>();
	double lhs_ = coerce_cast<double>(lhs);
	double rhs_ = coerce_cast<double>(rhs);
	if constexpr ( Method == Absolute )
		return lhs_ - rhs_;
	else if constexpr ( Method == RefLhs )
		return (lhs_ - rhs_) / lhs_;
	else if constexpr ( Method == RefRhs )
		return (lhs_ - rhs_) / rhs_;
	else
		static_assert(dependent_false<L>, "unsupported difference method");
}

// Compute signed absolute or relative difference
// - Use relative comparison if relative=true
// - For relative diff, referent determines the reference used
template<Num L, Num R>
double diff(
	const L query_v,
	const R table_v,
	const bool relative,
	const Ref referent = Query) noexcept
{
	if ( relative )
		switch(referent) {
			case Query: return diff<RefLhs>(query_v, table_v);
			case Table: return diff<RefRhs>(query_v, table_v);
		}
	else
		return diff<Absolute>(query_v, table_v);
}

// Does x neighbor ref within some tolerance(s)?
// - Tolerances and whether to use relative comparison are both per-dimension
// - For dimensions using relative diff, reference determines reference used
// - All dimensions using relative diff use the same referent
template<Num T = double, Vec L, Vec R, Vec Tol, Vec Rel>
bool near(
	const L query_v,
	const R table_v,
	const Tol tolerance,
	const Rel relative,
	const Ref referent = Query) noexcept
{
	assert(query_v.ssize() == table_v.ssize());
	assert(query_v.ssize() == tolerance.ssize());
	assert(query_v.ssize() == relative.ssize());
	for ( ptrdiff_t i = 0; i < query_v.ssize(); ++i )
	{
		double dx = diff(
			coerce_cast<T>(query_v[i]),
			coerce_cast<T>(table_v[i]),
			coerce_cast<bool>(relative[i]),
			referent);
		if ( std::fabs(dx) > tolerance[i] )
			return false;
	}
	return true;
}

//// Binary search
//-----------------

// Binary search for query in table
// - Values of table MUST be sorted (duplicates are ok)
// - Differences <= tolerance are considered matches
// - If relative == true, then referent determines the reference
// - Default nomatch chosen so nomatch << 0 for signed types
template<Num Index = ptrdiff_t, Num T, Vec V>
Index bsearch(
	const T query,
	const V table,
	const double tolerance = 0,
	const bool relative = false,
	const Ref referent = Query,
	const Index nomatch = na_value<Index>())
{
	if ( table.len == 0 )
		return nomatch;
	Index lo = 0;
	Index hi = table.len - 1;
	while ( lo <= hi )
	{
		Index mid = (lo + hi) / 2;
		double dx = diff(query, table[mid], relative, referent);
		if ( dx > 0 )
			lo = mid + 1;
		else if ( dx < 0 )
			hi = mid - 1;
		else
			return mid;
	}
	double dlo = std::fabs(diff(query, table[lo], relative, referent));
	double dhi = std::fabs(diff(query, table[hi], relative, referent));
	if ( dlo <= dhi && dlo <= tolerance )
		return lo;
	if ( dhi <= dlo && dhi <= tolerance )
		return hi;
	return nomatch;
}

// Binary search for multiple queries in ref
// - Values of table MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<Num Index = ptrdiff_t, Vec L, Vec R>
void bsearch(
	vec<Index> index,
	const L query,
	const R table,
	const double tolerance = 0,
	const bool relative = false,
	const Ref referent = Query,
	const Index nomatch = na_value<Index>())
{
	for ( ptrdiff_t i = 0; i < query.len; ++i )
	{
		if ( is_na(query[i]) )
			index[i] = nomatch;
		else
		{
			index[i] = bsearch(
				query[i], 
				table, 
				tolerance, 
				relative, 
				referent,
				nomatch);
		}
	}
}

// Ranged binary search for query in table
// - Values of table MUST be sorted (duplicates are ok)
// - Differences <= tolerance are considered matches
// - If relative == true, then referent determines the reference
// - Returns bounds [start, stop) range of matches
template<Num Index = ptrdiff_t, Num T, Vec V>
bounds rsearch(
	const T query,
	const V table,
	const double tolerance = 0,
	const bool relative = false,
	const Ref referent = Query,
	const Index nomatch = na_value<Index>())
{
	Index i = bsearch(query, table, tolerance, relative, referent, nomatch);
	if ( 0 <= i && i < table.ssize() )
	{
		Index lo = i;
		while ( 0 < lo && lo < table.ssize() - 1 && 
			std::fabs(diff(query, table[lo - 1], relative, referent))
				<= tolerance )
		{
			--lo;
		}
		Index hi = i;
		while ( 0 < hi && hi < table.ssize() - 1 &&
			std::fabs(diff(query, table[hi + 1], relative, referent))
				<= tolerance )
		{
			++hi;
		}
		return {lo, hi + 1};
	}
	else
		return {nomatch, nomatch};
}

// Ranged binary search for multiple queries in ref
// - Values of table MUST be sorted (duplicated are accepted)
// - Differences <= tolerance are considered matches
// - Default nomatch chosen so nomatch << 0 for signed types
template<Num Index = ptrdiff_t, Vec L, Vec R>
void rsearch(
	vec<Index> start,
	vec<Index> stop,
	const L query,
	const R table,
	const double tolerance = 0,
	const bool relative = false,
	const Ref referent = Query,
	const Index nomatch = na_value<Index>())
{
	for ( ptrdiff_t i = 0; i < query.len; ++i )
	{
		if ( is_na(query[i]) )
		{
			start[i] = nomatch;
			stop[i] = nomatch;
		}
		else
		{
			bounds b = rsearch(
				query[i], 
				table, 
				tolerance, 
				relative, 
				referent,
				nomatch);
			start[i] = b.start;
			stop[i] = b.stop;
		}
	}
}

//// K-D search
//--------------

// KDTree for nearest neighbor and range searches
// - Owner is responsible for managing memory
// - Table rows are observations and cols are features
// - Tree builds vectors of indices of children
// - NA or negative indices indicate leaf nodes
template<Num Index, Num T>
struct kdtree
{
	mat<T> table;
	vec<Index> left;
	vec<Index> right;
	Index root = na_value<Index>();

	ptrdiff_t ssize() const noexcept
	{
		return table.nrows();
	}

	bool has_left(Index node) const noexcept
	{
		return 0 <= left[node] && left[node] < table.nrows();
	}

	bool has_right(Index node) const noexcept
	{
		return 0 <= right[node] && right[node] < table.nrows();
	}

	// Build the tree and return the index of the root node
	Index build()
	{
		// invariants
		assert(left.len == table.nrows());
		assert(right.len == table.nrows());
		if ( table.ssize() <= 0 )
			return root;
		// initialize indices
		local_vec<Index> index{table.nrows()};
		index.fill_seq();
		left.fill_na();
		right.fill_na();
		// find root from median of first dim
		vec<T> column = table.col(0);
		qsort_index(index.borrow(), column);
		Index mid = table.nrows() / 2;
		// handle duplicates and update root
		while ( mid > 0 ) {
			if ( column.compare(index[mid - 1], index[mid]) == 0 )
				--mid;
			else
				break;
		}
		root = index[mid];
		// initialize stack
		Index top = -1;
		struct frame { Index parent, depth, start, stop; };
		auto stack = std::make_unique<frame[]>(max_depth(table.nrows()));
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
		if ( mid + 1 < table.nrows() ) {
			stack[++top] = {
				.parent = root,
				.depth = 1,
				.start = mid + 1,
				.stop = static_cast<Index>(table.nrows()),
			};
		}
		// recursively build the tree
		while ( top >= 0 )
		{
			// pop stack
			frame cur = stack[top--];
			// get current dim
			column = table.col(cur.depth % table.ncols());
			// find median of current dim within span of unprocessed rows
			if ( table.ncols() > 1 )
				qsort_index(index.borrow(), column, {cur.start, cur.stop});
			mid = (cur.start + cur.stop) / 2;
			while ( mid > cur.start ) {
				if ( column.compare(index[mid - 1], index[mid]) == 0 )
					--mid;
				else
					break;
			}
			// insert child under parent
			vec<T> previous = table.col((cur.depth - 1) % table.ncols());
			if ( previous.compare(index[mid], cur.parent) < 0 )
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

	// Apply Callable to indices of table within tolerance(s) of query
	// - Callable will be called as f(Index)
	// - Both tolerance and relative parameters are per-dimension
	// - Returns the count of hits
	template<UnaryOp Callable, Vec V, Vec Tol, Vec Rel>
	Index range_apply(
		Callable f,
		const V query,
		const Tol tolerance,
		const Rel relative,
		const Ref referent = Query) const
	{
		// invariants
		assert(query.ssize() == table.ncols());
		assert(tolerance.ssize() == table.ncols());
		assert(relative.ssize() == table.ncols());
		// initialize stack
		Index top = -1;
		struct frame { Index node, depth; };
		auto stack = std::make_unique<frame[]>(max_depth(table.nrows()));
		stack[++top] = { root, 0 };
		// initialize count of hits
		Index count = 0;
		// recursively search tree
		while ( top >= 0 )
		{
			// pop node
			frame cur = stack[top--];
			Index i = cur.depth % table.ncols();
			// compute distances
			double ds = diff(
				query[i], 
				table[{cur.node, i}], 
				relative[i], 
				referent);
			double du = std::fabs(ds);
			// is this a hit?
			bool is_hit = near(
				query, 
				table.row(cur.node), 
				tolerance, 
				relative, 
				referent);
			if ( is_hit )
			{
				f(cur.node);
				++count;
			}
			// search left subtree?
			if ( has_left(cur.node) && (ds < 0 || du <= tolerance[i]) )
				stack[++top] = { left[cur.node], cur.depth + 1 };
			// search right subtree?
			if ( has_right(cur.node) && (ds > 0 || du <= tolerance[i]) )
				stack[++top] = { right[cur.node], cur.depth + 1 };
		}
		return count;
	}

	// Get count of rows in table within tolerance of query
	template<Vec V, Vec Tol, Vec Rel>
	Index range_count(
		const V query,
		const Tol tolerance,
		const Rel relative,
		const Ref referent = Query) const
	{
		return range_apply(
			sink<Index>{},
			query,
			tolerance,
			relative,
			referent);
	}

	// Find indices of rows in table within tolerance of query
	// - Fills index with hits (up to index.len)
	// - Results ordering is arbitary
	// - Returns count of hits (may be > index.len)
	template<Vec V, Vec Tol, Vec Rel>
	Index range_search(
		vec<Index> index,
		const V query,
		const Tol tolerance,
		const Rel relative,
		const Ref referent = Query) const
	{
		return range_apply(
			sink<Index>{index},
			query,
			tolerance,
			relative,
			referent);
	}

	// Find indices of the K-nearest neighbors of a query in table
	// - Where K == index.len == dists.len
	// - Fills index with hits
	// - Fills dists with distances to hits
	// - Results ordered according to distance
	template<Vec V>
	void knn_search(
		vec<Index> index,
		vec<double> dists,
		const V query,
		const Norm p = L2) const
	{
		// invariants
		assert(query.ssize() == table.ncols());
		assert(index.ssize() == dists.ssize());
		// initialize stack
		Index top = -1;
		struct frame { Index node, depth; };
		auto stack = std::make_unique<frame[]>(max_depth(table.nrows()));
		stack[++top] = { root, 0 };
		// initialize KNN
		Index k = index.ssize();
		index.fill_na();
		dists.fill(huge_positive_value<double>());
		// recursively search tree
		while ( top >= 0 )
		{
			// pop node
			frame cur = stack[top--];
			Index i = cur.depth % table.ncols();
			// compute distances
			double ds = diff(query[i], table[{cur.node, i}]);
			double du = std::fabs(ds);
			double D = dist<double>(query, table.row(cur.node), p);
			// is this a hit?
			if ( D <= dists[k - 1] )
			{
				Index j = k - 1;
				// process strictly better nodes and/or ties
				if ( D <= dists[j] || cur.node < index[j] )
				{
					index[j] = cur.node;
					dists[j] = D;
					// sort neighbor into place
					while ( j > 0 && dists.compare(j, j - 1) <= 0 )
					{
						// break ties by index
						if ( dists.compare(j, j - 1) < 0 
							|| index.compare(j, j - 1) < 0 )
						{
							dists.swap(j, j - 1);
							index.swap(j, j - 1);
							--j;
						}
						else
							break;
					}
				}
			}
			// search left subtree?
			if ( has_left(cur.node) && (ds < 0 || du <= dists[k - 1]) )
				stack[++top] = { left[cur.node], cur.depth + 1 };
			// search right subtree?
			if ( has_right(cur.node) && (ds > 0 || du <= dists[k - 1]) )
				stack[++top] = { right[cur.node], cur.depth + 1 };
		}
	}

	#ifdef USING_R
	static kdtree<Index,T> from(SEXP obj)
	{
		return {
			.table = r_mat<T>(VECTOR_ELT(obj, 0)),
			.left = r_vec<Index>(VECTOR_ELT(obj, 1)),
			.right = r_vec<Index>(VECTOR_ELT(obj, 2)),
			.root = *data_ptr<Index>(VECTOR_ELT(obj, 3)),
		};
	}
	#endif // USING_R
};

// Range count kernel
template<Num Index, Num T, Vec Tol, Vec Rel>
struct range_counts
{
	vec<Index> counts;    // out
	kdtree<Index,T> tree; // in
	mat<T> query;         // in
	Tol tolerance;        // in
	Rel relative;         // in
	Ref referent;         // in

	ptrdiff_t ssize() const
	{
		return query.nrows();
	}

	void operator()(bounds b, task ctx)
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
		{
			counts[i] = tree.range_count(
				query.row(i),
				tolerance,
				relative,
				referent);
		}
	}
};

// Range search kernel
template<Num Index, Num T, Vec Tol, Vec Rel>
struct range_searches
{
	vecs_pack<Index,Index> index; // out
	kdtree<Index,T> tree;         // in
	mat<T> query;                 // in
	Tol tolerance;                // in
	Rel relative;                 // in
	Ref referent;                 // in

	ptrdiff_t ssize() const { return query.nrows(); }

	void operator()(bounds b, task ctx)
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
		{
			tree.range_search(
				index.get(i),
				query.row(i),
				tolerance,
				relative,
				referent);
			qsort(index.get(i));
		}
	}
};

// KNN search kernel
template<Num Index, Num T>
struct knn_searches
{
	vecs_pack<Index,Index> index;  // out
	vecs_pack<double,Index> dists; // out
	kdtree<Index,T> tree;          // in
	mat<T> query;                  // in
	Norm p;                        // in

	ptrdiff_t ssize() const { return query.nrows(); }

	void operator()(bounds b, task ctx)
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
		{
			tree.knn_search(
				index.get(i),
				dists.get(i),
				query.row(i),
				p);
		}
	}
};

#endif // CARDINAL_CORE_SEARCH
