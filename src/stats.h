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
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
		out_sums[i] = kern_sum<T>(x.col_vctr(i));
}

template<typename T>
void col_sums(
	const matrix<T> x, 
	double * out_sums,
	int num_threads = 1)
{
	if ( x.nrows == 0 || x.ncols == 0 )
		return;
	fill_buffer<double>(out_sums, x.ncols);
	num_threads = MIN2(num_threads, x.ncols);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		chunks c = chunks{0, x.ncols, num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_sums<T>, x, c.next(), out_sums
			};
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
void colrange_scatter_sums(
	const matrix<T> x, 
	const slice range,
	const int * group,
	const size_t ngroups,
	double * out_sums)
{
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
	{
		double * out_sums_i = out_sums + (i * ngroups);
		kern_scatter_sum(x.col_vctr(i), group, ngroups, out_sums_i);
	}
}

template<typename T>
void col_scatter_sums(
	const matrix<T> x, 
	const int * group,
	const size_t ngroups, 
	double * out_sums,
	int num_threads = 1)
{
	if ( x.nrows == 0 || x.ncols == 0 )
		return;
	fill_buffer<double>(out_sums, x.ncols);
	num_threads = MIN2(num_threads, x.ncols);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		chunks c = chunks{0, x.ncols, num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				colrange_scatter_sums<T>, 
				x, 
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
		colrange_scatter_sums<T>(
			x, 
			x.all_cols(), 
			group,
			ngroups,
			out_sums);
	}
}

#endif // CARDINAL_CORE_STATS
