#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include <atomic>
#include <thread>
#include "core.h"

//// Kernels
//-----------
// Distribute kernels to runners for computation

// A Kernel supports distributed computation on a range of indices
// - MUST be trivially copyable as a struct
// - MUST implement ssize() -> ptrdiff_t
// - MUST implement operator()(bounds b)
template<class F>
concept Kernel = 
	std::is_standard_layout_v<F> &&
	std::is_trivially_copyable_v<F> &&
	requires (std::remove_cvref_t<F>& f, ptrdiff_t i, bounds b)
	{
		{ f.ssize() } -> std::convertible_to<ptrdiff_t>;
		{ f(b) };
	};

// Chunk items for processing
struct chunker
{
	ptrdiff_t nchunks;
	ptrdiff_t nitems;

	bounds operator()(const ptrdiff_t i) const noexcept
	{
		// return early on special cases
		if ( i < 0 )
			return {0, 0};
		if ( i >= nchunks )
			return {nitems, nitems};
		// number of items per chunk
		ptrdiff_t chunksize = nitems / nchunks;
		// leftover items remaining
		ptrdiff_t remainder = nitems % nchunks;
		// distribute remainder across chunks
		if ( i < remainder ) {
			// put a remainder in this chunk
			if ( remainder > 0 )
				++chunksize;
			return {
				.start = chunksize * i,
				.stop = chunksize * (i + 1),
			};
		}
		else {
			// don't put a remainder in this chunk
			ptrdiff_t offset = (chunksize + 1) * remainder;
			return {
				.start = offset + (chunksize * (i - remainder)),
				.stop = offset + (chunksize * (i - remainder + 1)),
			};
		}
	}
};

// Executes a Kernel on exclusive index ranges
template<Kernel F>
struct runner
{
	F kernel;
	std::atomic<ptrdiff_t> * counter;
	chunker irange;
	
	void operator()()
	{
		do {
			ptrdiff_t i = counter->fetch_add(1, std::memory_order_relaxed);
			if ( 0 <= i && i < irange.nchunks )
				kernel(irange(i));
			else
				return;
		}
		while (true);
	}
};

// Dispatch a Kernel to parallel runners
// - Each thread's runner executes work chunked into ntasks
// - Runs on *this* thread if nthreads = 0
struct dispatcher
{
	std::atomic<ptrdiff_t> counter;
	std::thread * workers;
	int nthreads;
	int ntasks;
	bool active;

	dispatcher(int nthreads = 1, int ntasks = 1) : 
		workers(new std::thread[nthreads]),
		nthreads(nthreads),
		ntasks(ntasks),
		active(false) {}

	~dispatcher()
	{
		collect();
		delete[] workers;
	}

	dispatcher(const dispatcher&) = delete;
	dispatcher(dispatcher&&) = delete;
	dispatcher& operator=(const dispatcher&) = delete;
	dispatcher& operator=(dispatcher&&) = delete;

	template<Kernel F>
	void dispatch(F kernel)
	{
		if ( !active )
		{
			active = true;
			counter.store(0, std::memory_order_relaxed);
			chunker irange = { .nchunks = ntasks, .nitems = kernel.ssize() };
			if ( nthreads > 0 )
				for ( int i = 0; i < nthreads; ++i )
					workers[i] = std::thread{runner{kernel, &counter, irange}};
			else
				runner{kernel, &counter, irange}();
		}
	}

	void collect() noexcept
	{
		if ( active )
		{
			for ( int i = 0; i < nthreads; ++i )
				if ( workers[i].joinable() )
					workers[i].join();
			active = false;
		}
	}

	void stop() noexcept
	{
		counter.store(-1, std::memory_order_relaxed);
		collect();
	}
};

// Compute a Kernel over indices
// - Kernels receive mutually exclusive [start, stop) bounds
// - Kernels MUST exclusively modify data within these bounds
// - Limit nthreads and ntasks to the size of the kernel
// - If ntasks < nthreads, sets ntasks = nthreads
// - Runs on *this* thread if nthreads = 0
template<Kernel F>
void compute(F kernel, int nthreads = 1, int ntasks = 1)
{
	if ( nthreads < 0 )
		nthreads = 0;
	if ( nthreads > kernel.ssize() )
		nthreads = kernel.ssize();
	if ( ntasks < nthreads )
		ntasks = (nthreads == 0) ? 1 : nthreads;
	dispatcher work{nthreads, ntasks};
	work.dispatch(kernel);
}

//// Matrix statistics
//---------------------

template<typename T>
struct col_sums
{
	vec<double> sums;
	mat<T> x;

	ptrdiff_t ssize() const { return sums.len; }

	void operator()()
	{
		if ( x.row_stride > x.col_stride )
			for ( ptrdiff_t i = 0; i < x.nrows; ++i )
				sums += na_rm(x.row(i));
		else
			for ( ptrdiff_t j = 0; j < x.ncols; ++j )
				sums[j] = reduce<Add>(na_rm(x.col(j)));
	}

	void operator()(bounds b)
	{
		col_sums<T>{sums.slice(b), x.slice_cols(b)}();
	}
};

//// Test expressions
//--------------------

template<typename T>
void test_expression(vec<T> result, const vec<T> x, const vec<int> index)
{
	auto _x = na_rm(x);
	result.assign(gather(index, ufunc<Log1p>(_x + _x)));
}

#endif // CARDINAL_CORE_KERNELS
