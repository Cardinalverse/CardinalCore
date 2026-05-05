#ifndef CARDINAL_CORE_STAT
#define CARDINAL_CORE_STAT

#include <thread>
#include "prelude.h"

//// Unary kernel
//----------------

template<typename T>
double do_sum(
	const T * x,
	const size_t len,
	const ptrdiff_t stride = 1,
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	double sum = 0;
	for ( size_t i = 0; i < len; ++i )
	{
		double xi = static_cast<double>(x[i * stride]);
		if ( isIncomparable(xi) )
			continue;
		if ( weights != nullptr )
			xi *= weights[i];
		if ( abs )
			xi = std::fabs(xi);
		sum += std::pow(xi, p);
	}
	return sum;
}

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
void col_sums(
	const matrix<T> x, 
	const int num_threads, 
	double * out_sums)
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
			workers[i] = std::thread{colrange_sums<T>, x, begin, end, out_sums};
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

#endif // CARDINAL_CORE_STAT
