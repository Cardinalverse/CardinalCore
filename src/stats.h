#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include <thread>
#include "kernels.h"

//// Matrix statistics
//---------------------

template<typename T>
void colrange_sums(
	const matrix<T> x, 
	const size_t begin, 
	const size_t end, 
	double * out_sums)
{
	for ( size_t col = begin; col < end; ++col )
	{
		out_sums[col] = do_sum(x.col(col), x.nrows, x.row_stride);
	}
}

template<typename T>
void colrange_sums_grouped(
	const matrix<T> x, 
	const size_t begin, 
	const size_t end, 
	const int * group,
	const size_t ngroups,
	double * out_sums)
{
	for ( size_t col = begin; col < end; ++col )
	{
		double * out_sums_col = out_sums + (col * ngroups);
		do_sum_grouped(x.col(col), x.nrows,
			out_sums_col, group, ngroups, x.row_stride);
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
		size_t begin = 0, end = chunksize;
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_sums<T>, x, begin, end, out_sums
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
		colrange_sums<T>(x, 0, x.ncols, out_sums);
	}
}

template<typename T>
void col_sums_grouped(
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
		size_t begin = 0, end = chunksize;
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_sums_grouped<T>, 
				x, 
				begin, 
				end, 
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
		colrange_sums_grouped<T>(
			x, 
			0, 
			x.ncols, 
			group,
			ngroups,
			out_sums);
	}
}

#endif // CARDINAL_CORE_STATS
