#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include <atomic>
#include <thread>
#include "core.h"

//// Kernels
//-----------
// Distribute kernels to runners for computation

// A Kernel supports distributed computation
// - MUST be trivially copyable as a struct
// - MUST implement ssize() -> ptrdiff_t
// - MUST implement operator()(ptrdiff_t i)
template<class F>
concept Kernel = 
	std::is_standard_layout_v<F> &&
	std::is_trivially_copyable_v<F> &&
	requires (std::remove_cvref_t<F>& f, ptrdiff_t i)
	{
		{ f.ssize() } -> std::convertible_to<ptrdiff_t>;
		{ f(i) };
	};

// Execute a Kernel on an exclusive index
template<Kernel F>
struct runner
{
	F kernel;
	std::atomic<ptrdiff_t> * counter;
	
	void operator()()
	{
		do {
			ptrdiff_t i = counter->fetch_add(1, std::memory_order_relaxed);
			if ( 0 <= i && i < kernel.ssize() )
				kernel(i);
			else
				return;
		}
		while (true);
	}
};

// Dispatch a Kernel to parallel runners
struct dispatcher
{
	std::atomic<ptrdiff_t> counter = 0;
	std::thread * workers;
	int nthreads;

	explicit dispatcher(int n) : 
		workers(new std::thread[n]), nthreads(n) {}

	~dispatcher()
	{
		collect();
		delete[] workers;
	}

	dispatcher(const dispatcher&) = delete;
	dispatcher(dispatcher&&) noexcept = delete;
	dispatcher& operator=(const dispatcher&) = delete;
	dispatcher& operator=(dispatcher&&) noexcept = delete;

	template<Kernel F>
	void apply(F kernel)
	{
		if ( nthreads > 0 )
			for ( int i = 0; i < nthreads; ++i )
				workers[i] = std::thread{runner{kernel, &counter}};
		else
			runner{kernel, &counter}();
	}

	void collect() noexcept
	{
		for ( int i = 0; i < nthreads; ++i )
			if ( workers[i].joinable() )
				workers[i].join();
	}
};

// Chunk items for processing
struct chunker
{
	ptrdiff_t nchunks;
	ptrdiff_t nitems;

	bounds operator()(const ptrdiff_t i) const noexcept
	{
		if ( i < 0 )
			return {0, 0};
		else if ( i == 0 && (nchunks <= 1 || nitems <= nchunks) )
			return {0, nitems};
		else if ( i >= nchunks )
			return {nitems, nitems};
		ptrdiff_t chunksize = nitems / nchunks;
		ptrdiff_t remainder = nitems % nchunks;
		if ( i < remainder ) {
			if ( remainder > 0 )
				++chunksize;
			return {
				.start = chunksize * i,
				.stop = chunksize * (i + 1),
			};
		}
		else {
			ptrdiff_t offset = (chunksize + 1) * remainder;
			return {
				.start = offset + (chunksize * (i - remainder)),
				.stop = offset + (chunksize * (i - remainder + 1)),
			};
		}
	}
};

// Compute a Kernel over indices
template<Kernel F>
void compute(F kernel, int nthreads = 1)
{
	if ( nthreads < 0 )
		nthreads = kernel.ssize();
	dispatcher mc{nthreads};
	mc.apply(kernel);
}

//// Matrix statistics
//---------------------

template<typename T>
struct col_sums
{
	vec<double> sums;
	mat<T> x;
	int nchunks = 1;

	ptrdiff_t ssize() const { return nchunks; }

	void operator()()
	{
		if ( x.row_stride > x.col_stride )
			for ( ptrdiff_t i = 0; i < x.nrows; ++i )
				sums += x.row(i);
		else
			for ( ptrdiff_t j = 0; j < x.ncols; ++j )
				sums[j] = reduce<Add>(x.col(j));
	}

	void operator()(ptrdiff_t i)
	{
		bounds b = chunker{nchunks, sums.len}(i);
		col_sums<T>{sums.slice(b), x.slice_cols(b)}();
	}
};

//// Test expressions
//--------------------

template<typename T>
void test_expression(vec<T> result, const vec<T> x, const vec<int> index)
{
	result.assign(gather(index, transform<Log1p>(transform<Add>(x, x))));
	result.assign(gather(index, transform<Log1p>(x + x)));
}

#endif // CARDINAL_CORE_KERNELS
