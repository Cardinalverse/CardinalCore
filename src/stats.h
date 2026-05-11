#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include <thread>
#include "kernels.h"

//// Matrix statistics
//---------------------

template<typename T>
void sliceapply_sums(
	const matrix<T> x,
	const Axis axis,
	const slice range,
	double * out_sums)
{
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
	{
		switch(axis) {
			case Rows:
				out_sums[i] = kern_sum<T>(x.row_vctr(i));
				break;
			case Columns:
				out_sums[i] = kern_sum<T>(x.col_vctr(i));
				break;
		}
	}
}

template<typename T>
void apply_sums(
	const matrix<T> x, 
	const Axis axis,
	double * out_sums,
	int num_threads = 1)
{
	if ( x.nrows == 0 || x.ncols == 0 )
		return;
	size_t n = x.dim(axis);
	fill_buffer<double>(out_sums, n);
	num_threads = MIN2(num_threads, n);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		chunks c = chunks{0, n, num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				sliceapply_sums<T>, x, axis, c.next(), out_sums
			};
		}
		for ( int i = 0; i < num_threads; ++i )
			workers[i].join();
		delete [] workers;
	}
	else
	{
		sliceapply_sums<T>(x, axis, x.all_along(axis), out_sums);
	}
}

template<typename T>
void sliceapply_scatter_sums(
	const matrix<T> x, 
	const Axis axis,
	const slice range,
	const int * group,
	const size_t ngroups,
	double * out_sums)
{
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
	{
		double * out_sums_i = out_sums + (ngroups * i);
		switch(axis) {
			case Rows:
				kern_scatter_sum(x.row_vctr(i), group, ngroups, out_sums_i);
				break;
			case Columns:
				kern_scatter_sum(x.col_vctr(i), group, ngroups, out_sums_i);
				break;
		}
	}
}

template<typename T>
void apply_scatter_sums(
	const matrix<T> x, 
	const Axis axis,
	const int * group,
	const size_t ngroups, 
	double * out_sums,
	int num_threads = 1)
{
	if ( x.nrows == 0 || x.ncols == 0 )
		return;
	size_t n = x.dim(axis);
	fill_buffer<double>(out_sums, n * ngroups);
	num_threads = MIN2(num_threads, n);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		chunks c = chunks{0, n, num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				sliceapply_scatter_sums<T>, 
				x, 
				axis,
				c.next(), 
				group,
				ngroups,
				out_sums
			};
		}
		for ( int i = 0; i < num_threads; ++i )
			workers[i].join();
		delete [] workers;
	}
	else
	{
		sliceapply_scatter_sums<T>(
			x, 
			axis,
			x.all_along(axis), 
			group,
			ngroups,
			out_sums);
	}
}

#endif // CARDINAL_CORE_STATS
