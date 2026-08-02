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

template<Num T = ptrdiff_t>
struct counter
{
	std::atomic<T> count;
	std::atomic<T> limit;

	T next() noexcept
	{
		return count.fetch_add(1, std::memory_order_relaxed);
	}

	T max() noexcept
	{
		return limit.load(std::memory_order_relaxed);
	}

	void reset(T n) noexcept
	{
		count.store(0, std::memory_order_relaxed);
		limit.store(n, std::memory_order_relaxed);
	}

	void stop() noexcept
	{
		limit.store(0, std::memory_order_relaxed);
	}
};

// Execute a Kernel on an exclusive index
template<Kernel F>
struct runner
{
	F kernel;
	counter<> * queue;
	
	void operator()()
	{
		do
		{
			ptrdiff_t i = queue->next();
			if ( 0 <= i && i < queue->max() )
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
	counter<> queue;
	std::thread * workers;
	int nthreads;

	explicit dispatcher(int n) : 
		workers(new std::thread[n]), nthreads(n) {}

	~dispatcher()
	{
		collect(false);
		delete[] workers;
	}

	dispatcher(const dispatcher&) = delete;
	dispatcher(dispatcher&&) noexcept = delete;
	dispatcher& operator=(const dispatcher&) = delete;
	dispatcher& operator=(dispatcher&&) noexcept = delete;

	template<Kernel F>
	void apply(F kernel)
	{
		collect(true);
		queue.reset(kernel.ssize());
		for ( int i = 0; i < nthreads; ++i )
			workers[i] = std::thread{runner{kernel, &queue}};
	}

	void collect(bool force = false) noexcept
	{
		if ( force )
			queue.stop();
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
		if ( nitems <= nchunks ) {
			if ( i == 1 )
				return {0, nitems};
			else
				return {0, 0};
		}
		else if ( nchunks <= 1 )
			return {0, nitems};
		else if ( nchunks <= i )
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
	if ( nthreads >= 1 ) {
		dispatcher mc{nthreads};
		mc.apply(kernel);
	}
	else {
		for ( int i = 0; i < nthreads; ++i )
			kernel(i);
	}
}

//// Matrix statistics
//---------------------

template<typename T>
struct col_sums
{
	vec<double> sums;
	mat<T> x;
	int nchunks = 1;

	ptrdiff_t ssize() const { return sums.len; }

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
