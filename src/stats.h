#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include "parallel.h"
#include "kernels.h"

//// Matrix statistics
//---------------------

template<typename Kernel, typename T>
struct kern_col_sums
{
	vec<double> out;
	const mat<T> x;
	const Kernel kern = {};

	void operator()()
	{
		if ( x.row_stride > x.col_stride )
		{
			// row-major
			for ( ptrdiff_t i = 0; i < x.nrows; ++i )
				elementwise<Add>(out, x.row(i), kern);
		}
		else
		{
			// col-major
			for ( ptrdiff_t j = 0; j < x.ncols; ++j )
				out[j] = reduce<Add>(x.col(j), kern);
		}
	}
};

template<typename Kernel, typename T>
void col_sums(
	vec<double> out,
	const mat<T> x, 
	const Kernel kern = {},
	int num_threads = 1)
{
	num_threads = MIN2(num_threads, x.ncols);
	if ( num_threads > 1 )
	{
		bool ok = true;
		threads work{num_threads, &ok};
		chunks chunk{num_threads, x.ncols};
		for ( int i = 0; i < num_threads; ++i )
		{
			slice s = chunk(i);
			work.tasks[i] = std::thread{
				kern_col_sums{out.subset(s), x.subset_cols(s), kern}
			};
		}
		if ( !work.join_all() )
			Rf_error("one or more threads failed");
	}
	else
	{
		kern_col_sums{out, x, kern}();
	}
}

#endif // CARDINAL_CORE_STATS
