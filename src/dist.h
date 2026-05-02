#ifndef CARDINAL_CORE_DIST
#define CARDINAL_CORE_DIST

#include "compare.h"
#include "utils.h"

//// Mikowski distance 
//---------------------

// computes Mikowski distances between vectors x and y
// * dist = sum(w_i * |x_i - y_i|^p)^(1/p)
// returns: the distance
template<typename T>
double dist_mkw(
	const T * x,
	const T * y,
	size_t len,      // length of x and y
	double p = 2,
	int x_stride = 1,
	int y_stride = 1,
	const double * weights = nullptr)
{
	double D = 0;
	double w_i = 1;
	for ( size_t i = 0; i < len; ++i )
	{
		if ( weights != nullptr )
			w_i = weights[i];
		else
			w_i = 1;
		double d_i = qdiff(x[i * x_stride], y[i * y_stride]);
		D += std::pow(std::fabs(d_i), p);
	}
	return std::pow(D, 1 / p);
}

// computes Mikowski distances on rows of matrices of x and y
// * dist = sum(w_i * |x_i - y_i|^p)^(1/p)
// returns: the distance
template<typename T>
void do_dist_mkw(
	const T * x,
	const T * y,
	size_t len,      // length of x and y
	double p = 2,
	int x_stride = 1,
	int y_stride = 1,
	const double * weights = nullptr)
{
	double D = 0;
	double w_i = 1;
	for ( size_t i = 0; i < len; ++i )
	{
		if ( weights != nullptr )
			w_i = weights[i];
		else
			w_i = 1;
		double d_i = qdiff(x[i * x_stride], y[i * y_stride]);
		D += std::pow(std::fabs(d_i), p);
	}
	return std::pow(D, 1 / p);
}

