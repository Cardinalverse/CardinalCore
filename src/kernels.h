#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include <atomic>
#include <thread>
#include "core.h"

//// Kernels
//-----------
// Distribute kernels to runners for computation

// Task context
struct task { int task_id, thread_id; };

// A Kernel supports distributed computation on index ranges
// - MUST be trivially copyable as a struct
// - MUST implement ssize() -> ptrdiff_t
// - MUST implement operator()(bounds b, task ctx)
template<class F>
concept Kernel = 
	std::is_standard_layout_v<F> &&
	std::is_trivially_copyable_v<F> &&
	requires (std::remove_cvref_t<F>& f, bounds b, task ctx)
	{
		{ f.ssize() } -> std::convertible_to<ptrdiff_t>;
		{ f(b, ctx) };
	};

// Plan for parallel processing
struct plan
{
	int nthreads;
	int ntasks;

	// Get a plan ensuring nthreads <= ntasks <= nitems
	plan operator()(const ptrdiff_t nitems) const noexcept
	{
		plan p = {.nthreads = nthreads, .ntasks = ntasks};
		// Clamp nthreads to [0, nitems]
		p.nthreads = p.nthreads < 0 ? 0 : 
			(p.nthreads > nitems ? nitems : p.nthreads);
		// Clamp ntasks to [1, nthreads]
		p.ntasks = p.ntasks < 1 ? 1 :
			(p.ntasks < p.nthreads ? p.nthreads : p.ntasks);
		return p;
	}
};

// Chunk items for processing
// - Partition nitems into ntasks
// - A task is a half-open range of indices [start, stop)
struct chunker
{
	ptrdiff_t nitems;
	int ntasks;

	// Get index range for task i
	bounds operator()(const ptrdiff_t i) const noexcept
	{
		// return early on special cases
		if ( i < 0 )
			return {0, 0};
		if ( i >= ntasks )
			return {nitems, nitems};
		// number of items per task
		ptrdiff_t tasksize = nitems / ntasks;
		// leftover items remaining
		ptrdiff_t remainder = nitems % ntasks;
		// distribute remaining items across tasks
		if ( i < remainder ) {
			// give a remainder to this task
			if ( remainder > 0 )
				++tasksize;
			return {
				.start = tasksize * i,
				.stop = tasksize * (i + 1),
			};
		}
		else {
			// don't give a remainder to this task
			ptrdiff_t offset = (tasksize + 1) * remainder;
			return {
				.start = offset + (tasksize * (i - remainder)),
				.stop = offset + (tasksize * (i - remainder + 1)),
			};
		}
	}
};

// Executes a Kernel on mutually exclusive index ranges
// - Increment atomic counter to get a new task_id
// - Get task's index range from the chunker
// - Continue until counter is exhausts ntasks
template<Kernel F>
struct runner
{
	F kernel;
	std::atomic<int> * counter;
	chunker part;
	
	void operator()(const int thread_id)
	{
		while (counter != nullptr)
		{
			int task_id = counter->fetch_add(1, std::memory_order_relaxed);
			if ( 0 <= task_id && task_id < part.ntasks )
			{
				task ctx = {
					.task_id = task_id, 
					.thread_id = thread_id,
				};
				kernel(part(task_id), ctx);
			}
			else
				return;
		};
	}
};

// Dispatch a Kernel to parallel runners
// - Each thread's runner executes work chunked into ntasks
// - Runs on *this* thread if nthreads == 0
struct dispatcher
{
	std::atomic<int> counter;
	std::thread * workers;
	int nthreads;
	int ntasks;
	bool active;

	dispatcher(int nthreads = 1, int ntasks = 1) : 
		workers(new std::thread[nthreads]),
		nthreads(nthreads),
		ntasks(ntasks),
		active(false) {}
	
	dispatcher(plan p) : dispatcher(p.nthreads, p.ntasks) {}

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
			chunker part = {
				.nitems = kernel.ssize(),
				.ntasks = ntasks,
			};
			if ( nthreads > 0 )
			{
				for ( int j = 0; j < nthreads; ++j )
					workers[j] = std::thread{runner{kernel, &counter, part}, j};
			}
			else
				runner{kernel, &counter, part}(0);
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

// Compute a Kernel
// - Limits nthreads and ntasks to the ssize() of the kernel
// - Runs on *this* thread if nthreads <= 0
template<Kernel F>
void compute(F kernel, plan p)
{
	dispatcher work{p(kernel.ssize())};
	work.dispatch(kernel);
}

template<Kernel F>
void compute(F kernel, int nthreads = 1, int ntasks = 1)
{
	plan p = {.nthreads = nthreads, .ntasks = ntasks};
	compute(kernel, p);
}

//// Matrix statistics
//---------------------

template<Mat M>
struct col_sums
{
	vec<double> sums{};
	M x{};

	ptrdiff_t ssize() const { return sums.len; }

	void operator()()
	{
		if ( x.prefer_rows() )
			for ( ptrdiff_t i = 0; i < x.nrows(); ++i )
				sums += mask(x.row(i));
		else
			for ( ptrdiff_t j = 0; j < x.ncols(); ++j )
				sums[j] = sum(mask(x.col(j)));
	}

	void operator()(bounds b, task ctx)
	{
		col_sums<M>{sums.slice(b), x.slice_cols(b)}();
	}
};

//// Test expressions
//--------------------

template<typename T>
void test_expression(vec<T> result, const vec<T> x, const vec<int> index)
{
	result.assign(gather(index, log1p(mask(x) + mask(x))));
}

#endif // CARDINAL_CORE_KERNELS
