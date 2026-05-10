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

// computes sum(w_i * {x_i - m}^p)
// * if weights are null, then all w_i = 1
// * if abs = true, {} is absolute value
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
		if ( isIncomparable(x.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
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

// computes out_sums[g] = sum(w_i * {x_i - m_g}^p)
// where groups g are given by group[i]
// * if weights are null, then all w_i = 1
// * if abs = true, {} is absolute value
// * returns sums via out_sums
template<typename T>
void kern_scatter_sum(
	const vctr<T> x,
	const int * group,
	const size_t ngroups,
	double * out_sums,
	const double * m = nullptr,
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	fill_buffer<double>(out_sums, ngroups);
	for ( size_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
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

// computes sum(w_i * {x_i @ y_i}^p)
// where operation @ is a Binop given by <Op>
// * if weights are null, then all w_i = 1
// * if abs = true, {} is absolute value
// returns: the sum
template<typename Tx, typename Ty, int Op>
double kern2_sum(
	const vctr<Tx> x,
	const vctr<Ty> y,
	const double * weights = nullptr,
	const bool abs = false,
	const double p = 1)
{
	double sum = 0;
	for ( size_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) || isIncomparable(y.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
		double yi = static_cast<double>(y.at(i));
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
