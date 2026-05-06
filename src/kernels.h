#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include "prelude.h"

//// Elementary operations
//-------------------------

#define OP_ADD 1
#define OP_SUB 2
#define OP_MUL 3
#define OP_DIV 4

template<int Op>
double do_binop(double x, double y);

template<> inline
double do_binop<OP_ADD>(double x, double y) { return x + y; }

template<> inline
double do_binop<OP_SUB>(double x, double y) { return x - y; }

template<> inline
double do_binop<OP_MUL>(double x, double y) { return x * y; }

template<> inline
double do_binop<OP_DIV>(double x, double y) { return x / y; }

//// Unary kernels
//-----------------

// computes sum(x_i^p) by default
// * if weights are given, computes sum(w_i * x_i^p)
// * if abs is true, computes sum(w_i * |x_i|^p)
// returns: the sum
template<typename T>
double do_sum(
	const T * x,
	const size_t len, // length of x, weights
	const ptrdiff_t stride = 1, // stride of x
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
		if ( abs )
			xi = std::fabs(xi);
		if ( weights != nullptr )
			sum += weights[i] * std::pow(xi, p);
		else
			sum += std::pow(xi, p);
	}
	return sum;
}

// computes out_sums[g] = sum((x_ig)^p) by default
// where groups g are given by group[i]
// * if weights are given, computes sum(w_i * (x_ig)^p)
// * if abs is true, computes sum(w_i * |x_ig|^p)
// * returns sums via out_sums
template<typename T>
void do_sum_grouped(
	const T * x,
	const size_t len, // length of x, group, weights
	double * out_sums,
	const int * group,
	const size_t ngroups, // length of out_sums
	const ptrdiff_t stride = 1, // stride of x
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	fill_buffer<double>(out_sums, ngroups);
	for ( size_t i = 0; i < len; ++i )
	{
		double xi = static_cast<double>(x[i * stride]);
		if ( isIncomparable(xi) )
			continue;
		if ( abs )
			xi = std::fabs(xi);
		if ( weights != nullptr )
			out_sums[group[i]] += weights[i] * std::pow(xi, p);
		else
			out_sums[group[i]] += std::pow(xi, p);
	}
}

//// Binary kernels
//------------------

// computes sum((x_i o y_i)^p) by default
// where operation o is given by <Op> e.g., OP_SUB or OP_MUL
// * if weights are given, computes sum(w_i * (x_i o y_i)^p)
// * if abs is true, computes sum(w_i * |x_i o y_i|^p)
// returns: the sum
template<typename Tx, typename Ty, int Op>
double do_sum2(
	const Tx * x,
	const Ty * y,
	const size_t len, // length of x, y, weights
	const ptrdiff_t x_stride = 1, // stride of x
	const ptrdiff_t y_stride = 1, // stride of y
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	double sum = 0;
	for ( size_t i = 0; i < len; ++i )
	{
		double xi = static_cast<double>(x[i * x_stride]);
		double yi = static_cast<double>(y[i * y_stride]);
		if ( isIncomparable(xi) || isIncomparable(yi) )
			continue;
		double di = do_binop<Op>(xi, yi);
		if ( abs )
			di = std::fabs(di);
		if ( weights != nullptr )
			sum += weights[i] * std::pow(di, p);
		else
			sum += std::pow(di, p);
	}
	return sum;
}

// computes out_sums[g] = sum((x_ig o y_g)^p) by default
// where operation o is given by <Op> e.g., OP_SUB or OP_MUL
// * if weights are given, computes sum(w_i * (x_ig o y_g)^p)
// * if abs is true, computes sum(w_i * |x_ig o y_g|^p)
// * returns sums via out_sums
template<typename Tx, typename Ty, int Op>
void do_sum2_grouped(
	const Tx * x,
	const Ty * y,
	const size_t len, // length of x, group, weights
	double * out_sums,
	const int * group,
	const size_t ngroups, // length of y, out_sums
	const ptrdiff_t x_stride = 1, // stride of y
	const ptrdiff_t y_stride = 1, // stride of y
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	fill_buffer<double>(out_sums, ngroups);
	for ( size_t i = 0; i < len; ++i )
	{
		double xi = static_cast<double>(x[i * x_stride]);
		double yg = static_cast<double>(y[group[i] * y_stride]);
		if ( isIncomparable(xi) || isIncomparable(yg) )
			continue;
		double dgi = do_binop<Op>(xi, yg);
		if ( abs )
			dgi = std::fabs(dgi);
		if ( weights != nullptr )
			out_sums[group[i]] += weights[i] * std::pow(dgi, p);
		else
			out_sums[group[i]] += std::pow(dgi, p);
	}
}

#endif // CARDINAL_CORE_KERNELS
