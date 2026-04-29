#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

//// Comparison
//--------------

inline bool isIncomparable(int x)
{
	return x == NA_INTEGER;
}

inline bool isIncomparable(double x)
{
	return ISNA(x) || ISNAN(x);
}

// compute signed absolute or relative difference
// * safe to use with incomparables (NAs and NaNs)
// * incomparables sort last/highest (NA >> Inf)
// * incomparables sort equal to each other (NA == NA)
// returns: the difference
template<typename T>
double qdiff(T x, T ref, bool relative_diff = false)
{
	if ( isIncomparable(x) && isIncomparable(ref) )
		return 0.0;
	else if ( isIncomparable(x) )
		return R_PosInf; // NAs sort last so (x - ref) => +Inf
	else if ( isIncomparable(ref) )
		return R_NegInf; // NAs sort last so (x - ref) => -Inf
	else
	{
		if ( relative_diff )
			return static_cast<double>(x - ref) / ref;
		else
			return static_cast<double>(x - ref);
	}
}

// apply qdiff over an n-length vector
template<typename T>
void do_qdiff_impl(
	const T * x,
	const T * ref,
	double * out_result,
	size_t n,
	bool relative_diff)
{
	for ( size_t i = 0; i < n; ++i )
	{
		out_result[i] = qdiff<T>(x[i], ref[i], relative_diff);
	}
}

//// Quick select/sort
//---------------------

#define LESSER(x, y) (qdiff((x), (y)) < 0)
#define GREATER(x, y) (qdiff((x), (y)) > 0)
#define NOT_EQUAL(x, y) (qdiff((x), (y)) != 0)
#define SWAP(x, y, T) do { T swap = x; x = y; y = swap; } while (false)

// initialize a buffer with values
template<typename T>
void fill_buffer(
	T * buffer,      // buffer to fill
	size_t size,     // size of buffer
	T init = 0,      // value to use to initialize
	T increment = 0, // increment init for each item
	int stride = 1)
{
	for ( size_t i = 0; i < size; i++ )
	{
		buffer[stride * i] = init;
		init += increment;
	}
}

// select a pivot and partition x around the pivot such that
// * partitions the _indices_ of x via out_index
// * all items left of pivot are <= pivot
// * all items right of pivot are >= pivot
// * incomparables sort last/highest (NA >> Inf)
// returns: pivot index
template<typename T, typename Index>
Index partition(
	const T * x,
	Index lo, // index of first item to consider in x
	Index hi, // index of last item to consider in x
	Index * out_index,
	bool init_out_index = false)
{
	// we get item k via x[at[k]]
	Index * at = out_index;
	// fill out_index with sequential indices
	if ( init_out_index )
		fill_buffer<Index>(out_index, hi - lo + 1, 0, 1);
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

// find the k-th ranked item of an array x
// * partially sorts the _indices_ of x via out_index
// * incomparables rank last/highest (NA >> Inf)
// returns: index of k-th item
template<typename T, typename Index, typename Rank>
T quick_select(
	T const * x,
	Index begin, // index of first item to consider
	Index end,   // one-past-the-end index
	Rank k,
	Index * out_index,
	bool init_out_index = false)
{
	// we get item k via x[at[k]]
	Index * at = out_index;
	// fill out_index with sequential indices
	if ( init_out_index )
		fill_buffer<Index>(out_index, end - begin, 0, 1);
	// recursively partition the array
	Index lo = begin, hi = end - 1;
	do {
		if ( lo == hi )
			return x[at[lo]];
		Index pivot = partition<T,Index>(x, lo, hi, at);
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

// find the k-th ranked items of an array for multiple k's
template<typename T, typename Index, typename Rank>
void do_qselect_impl(
	const T * x, 
	Index begin, 
	Index end, 
	const Rank * ks, 
	size_t ks_len,
	T * out_values)
{
	// validate input size
	size_t n;
	if ( end - begin > 0 )
		n = static_cast<size_t>(end - begin);
	else
		return;
	// set up working index buffer
	Index * work_index = R_Calloc(n, Index);
	fill_buffer<Index>(work_index, n, 0, 1);
	// loop through ks	
	for ( size_t i = 0; i < ks_len; ++i )
	{
		if ( i == 0 )
			out_values[0] = quick_select<T,Index,Rank>(x, 0, n, ks[0], work_index);
		else if ( ks[i] > ks[i - 1] )
			out_values[i] = quick_select<T,Index,Rank>(x, ks[i - 1] + 1, n, ks[i], work_index);
		else if ( ks[i] < ks[i - 1] )
			out_values[i] = quick_select<T,Index,Rank>(x, 0, ks[i - 1], ks[i], work_index);
		else
			out_values[i] = out_values[i - 1];
	}
	R_Free(work_index);
}

// sort an array x using Hoare's quicksort algorithm
// * sorts the _indices_ of x via out_index
// * incomparables rank last/highest (NA >> Inf)
template<typename T, typename Index>
void quick_sort(
	const T * x, 
	Index begin, 
	Index end, 
	Index * out_index,
	bool init_out_index = false,
	int linear_threshold = 8)
{
	// validate input size
	size_t n;
	if ( end - begin > 0 )
		n = static_cast<size_t>(end - begin);
	else
		return;
	// fill out_index with sequential indices
	if ( init_out_index )
		fill_buffer<Index>(out_index, n, 0, 1);
	// initialize the stack
	size_t stack_size = 2 * std::ceil(std::log2(n) + 1);
	Index * stack = R_Calloc(stack_size, Index);
	Index top = -1;
	Index lo = begin, hi = end - 1;
	stack[++top] = lo;
	stack[++top] = hi;
	// we get item k via x[at[k]]
	Index * at = out_index;
	// recursively partition the array
	while ( top >= 0 )
	{
		// pop and partition current subarray
		hi = stack[top--];
		lo = stack[top--];
		if ( hi - lo < linear_threshold )
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
		Index pivot = partition<T,Index>(x, lo, hi, at);
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

#endif // CARDINAL_CORE_SEARCH
