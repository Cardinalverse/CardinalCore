#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include "parallel.h"
#include "kernels.h"

//// Matrix statistics
//---------------------

template<typename T, int Tform, int Op>
void apply_kern(
	const matrix<T> x,
	const Axis axis,
	const slice range,
	double * out_values)
{
	for ( ptrdiff_t i = range.begin; i < range.end; ++i )
	{
		switch(axis) {
			case Rows:
				out_values[i] = kern_reduce<T,Tform,Op>(x.row_vctr(i));
				break;
			case Columns:
				out_values[i] = kern_reduce<T,Tform,Op>(x.col_vctr(i));
				break;
		}
	}
}

template<typename T, int Tform, int Op>
void par_apply_kern(
	const matrix<T> x, 
	const Axis axis,
	double * out_values,
	int num_threads = 1)
{
	num_threads = MIN2(num_threads, x.dim(axis));
	if ( num_threads > 1 )
	{
		bool ok = true;
		threads work{num_threads, &ok};
		chunks c{0, x.dim(axis), num_threads};
		for ( int i = 0; i < num_threads; ++i )
		{
			work.tasks[i] = std::thread{
				apply_kern<T,Tform,Op>, 
				x, axis, c.next(), out_values
			};
		}
		if ( !ok )
			Rf_error("one or more threads failed");
	}
	else
	{
		apply_kern<T,Tform,Op>(
			x, axis, x.all_along(axis), out_values);
	}
}

#endif // CARDINAL_CORE_STATS
