#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include "prelude.h"

//// Elementary operations
//-------------------------

enum Binop {
	Addition,
	Subtraction,
	Multiplication,
	Division
};

template<int Op>
double do_binop(double x, double y);

template<> inline
double do_binop<Addition>(double x, double y) { return x + y; }

template<> inline
double do_binop<Subtraction>(double x, double y) { return x - y; }

template<> inline
double do_binop<Multiplication>(double x, double y) { return x * y; }

template<> inline
double do_binop<Division>(double x, double y) { return x / y; }

//// Unary kernels
//-----------------

// computes sum((x_i - m)^p) by default
// * if weights are given, computes sum(w_i * (x_i - m)^p)
// * if abs is true, computes sum(w_i * |x_i - m|^p)
// returns: the sum
template<typename T>
double kern_sum(
	const vctr<T> x,
	const double m = 0,
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	double sum = 0;
	for ( size_t i = 0; i < x.len; ++i )
	{
		double xi = static_cast<double>(x.at(i));
		if ( isIncomparable(xi) )
			continue;
		if ( m != 0 )
			xi -= m;
		if ( abs )
			xi = std::fabs(xi);
		if ( weights != nullptr )
			sum += weights[i] * std::pow(xi, p);
		else
			sum += std::pow(xi, p);
	}
	return sum;
}

// computes out_sums[g] = sum((x_i - m_g)^p) by default
// where groups g are given by group[i]
// * if weights are given, computes sum(w_i * (x_i - m_g)^p)
// * if abs is true, computes sum(w_i * |x_i - m_g|^p)
// * returns sums via out_sums
template<typename T>
void kern_sum_grouped(
	const vctr<T> x,
	double * out_sums,
	const int * group,
	const size_t ngroups, // length of out_sums
	const double * m = nullptr,
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	fill_buffer<double>(out_sums, ngroups);
	for ( size_t i = 0; i < x.len; ++i )
	{
		double xi = static_cast<double>(x.at(i));
		if ( isIncomparable(xi) )
			continue;
		if ( m != nullptr )
			xi -= m[group[i]];
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
double kern2_sum(
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

#endif // CARDINAL_CORE_KERNELS
