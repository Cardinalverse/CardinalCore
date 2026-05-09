#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include <thread>
#include "kernels.h"

//// Matrix statistics
//---------------------

template<typename T>
void colrange_sums(
	const matrix<T> x, 
	const slice range,
	double * out_sums)
{
	for ( size_t i = range.begin; i < range.end; ++i )
		out_sums[i] = kern_sum<T>(x.col_vctr(i));
}

template<typename T>
void colrange_scatter_sums(
	const matrix<T> x, 
	const slice range,
	const int * group,
	const size_t ngroups,
	double * out_sums)
{
	for ( size_t i = range.begin; i < range.end; ++i )
	{
		double * out_sums_i = out_sums + (i * ngroups);
		kern_sum_grouped(x.col_vctr(i), group, ngroups, out_sums_i);
	}
}

template<typename T>
void col_sums(
	const matrix<T> x, 
	double * out_sums,
	const int num_threads = 1)
{
	fill_buffer<double>(out_sums, x.ncols);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		int chunksize = x.ncols / num_threads;
		if ( x.ncols % num_threads > 0 )
			chunksize += 1;
		ptrdiff_t begin = 0, end = chunksize;
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_sums<T>, x, slice{begin, end}, out_sums
			};
			begin += chunksize;
			end += chunksize;
			if ( end > x.ncols )
				end = x.ncols;
		}
		for ( int i = 0; i < num_threads; ++i )
			workers[i].join();
		delete [] workers;
	}
	else
	{
		colrange_sums<T>(x, x.all_cols(), out_sums);
	}
}

template<typename T>
void col_scatter_sums(
	const matrix<T> x, 
	const int * group,
	const size_t ngroups, 
	double * out_sums,
	const int num_threads = 1)
{
	fill_buffer<double>(out_sums, x.ncols);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		int chunksize = x.ncols / num_threads;
		if ( x.ncols % num_threads > 0 )
			chunksize += 1;
		ptrdiff_t begin = 0, end = chunksize;
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_scatter_sums<T>, 
				x, 
				slice{begin, end}, 
				group,
				ngroups,
				out_sums
			};
			begin += chunksize;
			end += chunksize;
			if ( end > x.ncols )
				end = x.ncols;
		}
		for ( int i = 0; i < num_threads; ++i )
			workers[i].join();
		delete [] workers;
	}
	else
	{
		colrange_scatter_sums<T>(
			x, 
			x.all_cols(), 
			group,
			ngroups,
			out_sums);
	}
}

#endif // CARDINAL_CORE_STATS
