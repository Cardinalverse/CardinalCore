#ifndef CARDINAL_CORE_APPLY
#define CARDINAL_CORE_APPLY

#include <thread>

//// Parallel apply
//-------------------
// Apply a callable over chunks in parallel

// Array of threads
// - Threads are joined at the end of the scope
struct threads
{
	std::thread * tasks;
	int nthreads;

	explicit threads(int n) : 
		tasks(new std::thread[n]), nthreads(n) {}

	~threads()
	{
		collect();
		delete[] tasks;
	}

	threads(const threads&) = delete;
	threads(threads&&) noexcept = delete;
	threads& operator=(const threads&) = delete;
	threads& operator=(threads&&) noexcept = delete;

	bool collect() noexcept
	{
		bool ok = true;
		for ( int i = 0; i < nthreads; ++i )
		{
			if ( !tasks[i].joinable() )
				continue;
			try {
				tasks[i].join();
			} catch (...) {
				ok = false;
			}
		}
		return ok;
	}
};

// Chunk items for processing
// - The nchunks is the number of chunks and MUST be >= 0
// - The size is the number of items and MUST be >= 0
// - Callable yields [start, stop) for ith chunk
template<typename Index = ptrdiff_t>
struct chunks
{
	int nchunks;
	Index size;

	bounds operator()(const Index i) const noexcept
	{
		if ( i < 0 )
			return {0, 0};
		if ( i >= nchunks )
			return {size, size};
		if ( nchunks <= 1 || size <= nchunks )
			return {0, size};
		Index chunksize = size / nchunks;
		Index remainder = size % nchunks;
		if ( i < remainder ) {
			if ( remainder > 0 )
				++chunksize;
			return {
				.start = chunksize * i,
				.stop = chunksize * (i + 1),
			};
		}
		else {
			Index offset = (chunksize + 1) * remainder;
			return {
				.start = offset + (chunksize * (i - remainder)),
				.stop = offset + (chunksize * (i - remainder + 1)),
			};
		}
	}
};

//// Apply a callable Kernel to chunks in parallel
// - Each chunk gets its own thread
// - Each chunk gets mutually exclusive [start, stop) bounds
// - The Kernel must implement operator(bounds)
template<typename Kernel, typename Index = ptrdiff_t>
void chunk_apply(
	Kernel kernel,
	const int nchunks,
	const Index size,
	const bool parallel = true)
{
	if ( nchunks > 1 )
	{
		chunks chunk{nchunks, size};
		if ( parallel )
		{
			threads work{nchunks};
			for ( int i = 0; i < nchunks; ++i )
				work.tasks[i] = std::thread{kernel, chunk(i)};
		}
		else
		{
			for ( int i = 0; i < nchunks; ++i )
				kernel(chunk(i));
		}
	}
	else
	{
		kernel(bounds{0, size});
	}
}

////// Matrix statistics
//---------------------

template<typename T>
struct kern_col_sums
{
	vec<double> out;
	const mat<T> x;

	void operator()()
	{
		kern_unop<Identity> kernel{};
		if ( x.row_stride > x.col_stride )
			for ( ptrdiff_t i = 0; i < x.nrows; ++i )
				elementwise<Add>(out, x.row(i), kernel);
		else
			for ( ptrdiff_t j = 0; j < x.ncols; ++j )
				out[j] = reduce<Add>(x.col(j), kernel);
	}

	void operator()(bounds b)
	{
		kern_col_sums<T>{out.slice(b), x.slice_cols(b)}();
	}
};

#endif // CARDINAL_CORE_APPLY
