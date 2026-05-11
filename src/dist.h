#ifndef CARDINAL_CORE_DIST
#define CARDINAL_CORE_DIST

#include <thread>
#include "kernels.h"

//// Mikowski distance 
//---------------------

template<typename T>
void sliceapply_dists_mkw(
	const matrix<T> x,
	const matrix<T> y,
	const Axis axis,
	const slice range,
	double * out_dists)
{
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
	{
		double * out_dists_i = out_dists + (y.dim(axis) * i);
		for ( size_t j = 0; j < y.dim(axis); ++j )
		{
			switch(axis) {
				case Rows:
					out_dists_i[j] = kern2_sum<T,T,Minus>(
						x.row_vctr(i), y.row_vctr(j), nullptr, true, 2);
					break;
				case Columns:
					out_dists_i[j] = kern2_sum<T,T,Minus>(
						x.col_vctr(i), y.col_vctr(j), nullptr, true, 2);
					break;
			}
		}
	}
}

template<typename T>
void apply_dists_mkw(
	const matrix<T> x, 
	const matrix<T> y, 
	const Axis axis,
	double * out_dists,
	int num_threads = 1)
{
	if ( x.nrows == 0 || x.ncols == 0 )
		return;
	size_t n = x.dim(axis);
	fill_buffer<double>(out_dists, n);
	num_threads = MIN2(num_threads, n);
	if ( num_threads > 1 )
	{
		std::thread * workers = new std::thread[num_threads];
		chunks c = chunks{0, n, num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			workers[i] = std::thread{
				sliceapply_dists_mkw<T>, x, y, axis, c.next(), out_dists
			};
		}
		for ( int i = 0; i < num_threads; ++i )
			workers[i].join();
		delete [] workers;
	}
	else
	{
		sliceapply_dists_mkw<T>(x, y, axis, x.all_along(axis), out_dists);
	}
}

#endif // CARDINAL_CORE_DIST
