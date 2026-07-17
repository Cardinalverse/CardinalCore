#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include "prelude.h"

//// Safety
//----------
// All these functions should be considered unsafe
// The caller takes responsiblity for ensuring:
// - pointers are allocated and initialized
// - arrays lengths are appropriate
// - memory is later freed
// These functions do not allocate

//// Unary operations
//--------------------

enum Unop {
	Noop,
	Abs,
	Log,
	Log2,
	Log1p,
	Exp,
	Exp2,
	Expm1
};

template<int Op>
double do_unop(double x);

template<> inline
double do_unop<Noop>(double x) { return x; }

template<> inline
double do_unop<Abs>(double x) { return std::abs(x); }

template<> inline
double do_unop<Log>(double x) { return std::log(x); }

template<> inline
double do_unop<Log2>(double x) { return std::log2(x); }

template<> inline
double do_unop<Log1p>(double x) { return std::log1p(x); }

template<> inline
double do_unop<Exp>(double x) { return std::exp(x); }

template<> inline
double do_unop<Exp2>(double x) { return std::exp2(x); }

template<> inline
double do_unop<Expm1>(double x) { return std::expm1(x); }

//// Binary operations
//---------------------

enum Binop {
	Add,
	Subtract,
	Multiply,
	Divide,
	Max,
	Min
};

template<int Op>
double do_binop(double x, double y);

template<> inline
double do_binop<Add>(double x, double y) { return x + y; }

template<> inline
double do_binop<Subtract>(double x, double y) { return x - y; }

template<> inline
double do_binop<Multiply>(double x, double y) { return x * y; }

template<> inline
double do_binop<Divide>(double x, double y) { return x / y; }

template<> inline
double do_binop<Max>(double x, double y) { return x > y ? x : y; }

template<> inline
double do_binop<Min>(double x, double y) { return x < y ? x : y; }

template<int Op>
double init_binop_accum()
{
	switch(Op) {
		case Add:
		case Subtract:
			return 0;
		case Multiply:
		case Divide:
			return 1;
		case Max:
			return NEG_INF;
		case Min:
			return POS_INF;
	}
}

//// Unary kernels
//-----------------

// kern_reduce<T,Tform,Reduce> -> double
// => Reduce(..., w[i] * Tform(x[i] - c)^p)
// - Tform is a Unop
// - Reduce is a Binop
// - if weights are null, then all w[i] = 1
// - incomparables are skipped
// returns: scalar double result from the reduction
// e.g., if Reduce is Add, then returns the sum
template<typename T, int Tform, int Reduce>
double kern_reduce(
	const vctr<T> x,
	const double c = 0,
	const double p = 1,
	const double * weights = nullptr)
{
	double accum = init_binop_accum<Reduce>();
	for ( isize i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
		xi = std::pow(do_unop<Tform>(xi - c), p);
		if ( weights != nullptr )
			xi = weights[i] * xi;
		accum = do_binop<Reduce>(accum, xi);
	}
	return accum;
}

// kern_accum<T,Tform,Reduce> -> void
// => out_accum[i] = Reduce(out_accum[i], w[i] * Tform(x[i] - c)^p))
// - Tform is a Unop
// - Reduce is a Binop
// - out_accum MUST be initialized already
// - out_accum MUST have length equal to x.len
// - if weights are null, then all w[i] = 1
// - incomparables are skipped
// returns: nothing
template<typename T, int Tform, int Reduce>
void kern_accum(
	const vctr<T> x,
	double * out_accum,
	const double c = 0,
	const double p = 1,
	const double * weights = nullptr)
{
	for ( isize i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
		xi = std::pow(do_unop<Tform>(xi - c), p);
		if ( weights != nullptr )
			xi = weights[i] * xi;
		out_accum[i] = do_binop<Reduce>(out_accum[i], xi);
	}
}

// kern_scatter<T,Tform,Reduce> -> void
// => out_accum[g] = Reduce(Reduce, w[i] * Tform(x[i] - c[g])^p)
// where groups g are given by group.index[i]
// - Tform is a Unop
// - Reduce is a Binop
// - out_accum MUST be initialized already
// - out_accum MUST have size >= x.len
// - if c is not null, c MUST have size >= groups.ngroups
// - if weights are null, then all w[i] = 1
// - incomparables are skipped
// returns: nothing
template<typename T, int Tform, int Reduce>
void kern_scatter(
	const vctr<T> x,
	const groups group,
	double * out_accum,
	const double * c = nullptr,
	const double p = 1,
	const double * weights = nullptr)
{
	for ( isize i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x.at(i)) )
			continue;
		double xi = static_cast<double>(x.at(i));
		isize g = group.index[i];
		if ( c != nullptr )
			xi = std::pow(do_unop<Tform>(xi - c[g]), p);
		else
			xi = std::pow(do_unop<Tform>(xi), p);
		if ( weights != nullptr )
			xi = weights[i] * xi;
		out_accum[g] = do_binop<Reduce>(out_accum[i], xi);
	}
}

#endif // CARDINAL_CORE_KERNELS
