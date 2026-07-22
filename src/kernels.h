#ifndef CARDINAL_CORE_KERNELS
#define CARDINAL_CORE_KERNELS

#include <cmath>
#include "core.h"

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

template<int Op>
struct kern_unop
{
	double operator()(double x) const 
	{
		return do_unop<Op>(x);
	}
};

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
double init_accum()
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

//// Vector operations
//---------------------

template<int Reduce, typename Kernel, typename T>
void elementwise(
	vec<double> out,
	const vec<T> x,
	const Kernel kern = {})
{
	for ( ptrdiff_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x[i]) )
			continue;
		double xi = static_cast<double>(x[i]);
		out[i] = do_binop<Reduce>(out[i], kern(xi));
	}
}

template<int Reduce, typename Kernel, typename T>
double reduce(
	const vec<T> x,
	const Kernel kern = {})
{
	double out = init_accum<Reduce>();
	for ( ptrdiff_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x[i]) )
			continue;
		double xi = static_cast<double>(x[i]);
		out = do_binop<Reduce>(out, kern(xi));
	}
	return out;
}

template<int Reduce, typename Kernel, typename Index, typename T>
void scatter(
	vec<double> out,
	const Index * index,
	const vec<T> x,
	const Kernel kern = {})
{
	for ( ptrdiff_t i = 0; i < x.len; ++i )
	{
		if ( isIncomparable(x[i]) )
			continue;
		double xi = static_cast<double>(x[i]);
		out[index[i]] = do_binop<Reduce>(out[index[i]], kern(xi));
	}
}

#endif // CARDINAL_CORE_KERNELS
