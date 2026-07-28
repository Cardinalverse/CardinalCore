#ifndef CARDINAL_CORE_CORE
#define CARDINAL_CORE_CORE

#include <cassert>
#include <cstddef>
#include <cmath>
#include <limits>
#include <type_traits>

//// Utility
//-----------
// Common idioms

#define MIN2(x, y) ((x) < (y) ? (x) : (y))
#define MIN3(x, y, z) (MIN2(MIN2((x), (y)), (z)))
#define MAX2(x, y) ((x) > (y) ? (x) : (y))
#define MAX3(x, y, z) (MAX2(MAX2((x), (y)), (z)))

//// Data types
//--------------
// Coercion and incomparables

template<typename T>
constexpr T make_incomparable()
{
	if constexpr ( std::is_floating_point_v<T> )
		return std::numeric_limits<T>::quiet_NaN();
	else
		return std::numeric_limits<T>::lowest();
}

#ifdef USING_R
template<>
inline int make_incomparable<int>() { return NA_INTEGER; }
template<>
inline double make_incomparable<double>() { return NA_REAL; }
#endif // USING_R

template<typename T>
constexpr bool is_incomparable(T v)
{
	if constexpr ( std::is_floating_point_v<T> )
		return std::isnan(v);
	else
		return v == make_incomparable<T>();
}

#ifdef USING_R
template<>
inline bool is_incomparable<int>(int v) { return v == NA_INTEGER; }
template<>
inline bool is_incomparable<double>(double v) { return ISNAN(v); }
#endif // USING_R

template<typename Out, typename In>
Out coerce_cast(In v)
{
	if ( is_incomparable(v) )
		return make_incomparable<Out>();
	else
		return static_cast<Out>(v);
}

//// Data Pointer
//----------------
// Get a mutable data pointer from a managed object

#ifdef USING_R
template<typename T>
T * data_ptr(SEXP x);
template<> inline
int * data_ptr<int>(SEXP x) { return INTEGER(x); }
template<> inline
double * data_ptr<double>(SEXP x) { return REAL(x); }
#endif // USING_R

//// Unary operations
//--------------------
// Generic unary functors

enum Unop {
	Identity,
	Abs,
	Log,
	Log2,
	Log1p,
	Exp,
	Exp2,
	Expm1
};

template<Unop Op, typename T = double>
struct unop {
	T invoke(T x);
	T operator()(T x)
	{
		if ( is_incomparable(x) )
			return make_incomparable<T>();
		else
			return unop<Op,T>::invoke(x);
	}
};
template<typename T>
struct unop<Identity,T> { 
	T invoke(T x) { return x; }
};
template<typename T>
struct unop<Abs,T> {
	T invoke(T x) { return std::abs(x); } 
};
template<typename T>
struct unop<Log,T> {
	T invoke(T x) { return std::log(x); } 
};
template<typename T>
struct unop<Log2,T> {
	T invoke(T x) { return std::log2(x); }
};
template<typename T>
struct unop<Log1p,T> {
	T invoke(T x) { return std::log1p(x); }
};
template<typename T>
struct unop<Exp,T> {
	T invoke(T x) { return std::exp(x); } 
};
template<typename T>
struct unop<Exp2,T> {
	T invoke(T x) { return std::exp2(x); } 
};
template<typename T>
struct unop<Expm1,T> {
	T invoke(T x) { return std::expm1(x); } 
};

//// Binary operations
//---------------------
// Generic binary functors

enum Binop {
	LHS,
	RHS,
	Add,
	Subtract,
	Multiply,
	Divide,
	Max,
	Min
};

template<Binop Op, typename T = double>
struct binop {
	T invoke(T lhs, T rhs);
	T operator()(T lhs, T rhs)
	{
		if ( is_incomparable(lhs) || is_incomparable(rhs) )
			return make_incomparable<T>();
		else
			return binop<Op,T>::invoke(lhs, rhs);
	}
};
template<typename T>
struct binop<LHS,T> {
	T invoke(T lhs, T rhs) { return lhs; }
};
template<typename T>
struct binop<RHS,T> {
	T invoke(T lhs, T rhs) { return rhs; }
};
template<typename T>
struct binop<Add,T> {
	T invoke(T lhs, T rhs) { return lhs + rhs; }
};
template<typename T>
struct binop<Subtract,T> {
	T invoke(T lhs, T rhs) { return lhs - rhs; }
};
template<typename T>
struct binop<Multiply,T> {
	T invoke(T lhs, T rhs) { return lhs * rhs; }
};
template<typename T>
struct binop<Divide,T> {
	T invoke(T lhs, T rhs) { return lhs / rhs; }
};
template<typename T>
struct binop<Max,T> {
	T invoke(T lhs, T rhs) { return lhs > rhs ? lhs : rhs; }
};
template<typename T>
struct binop<Min,T> {
	T invoke(T lhs, T rhs) { return lhs < rhs ? lhs : rhs; }
};

template<Binop Op, typename T = double>
T binop_init()
{
	switch(Op) {
		case LHS:
		case RHS:
		case Add:
		case Subtract:
			return coerce_cast<T>(0);
		case Multiply:
		case Divide:
			return coerce_cast<T>(1);
		case Max:
			return coerce_cast<T>(-INFINITY);
		case Min:
			return coerce_cast<T>(+INFINITY);
	}
}

//// Kernel operations
// -------------------
// Apply functions and kernels

template<typename Vec, typename Reduce, typename T>
T reduce(Vec input, Reduce op, T init) 
{
	T output = init;
	for ( size_t i = 0; i < input.size(); ++i )
		output = op(output, coerce_cast<T>(input[i]));
	return output;
}

template<Binop Op, typename Vec, typename T>
T reduce(Vec input) 
{
	return reduce(input, binop<Op,T>{}, init_binop<Op,T>());
}

//// Structs
//-----------
// Containers for common object types

// Index bounds
// - The interval is half-open: [start, stop)
struct bounds 
{
	ptrdiff_t start;
	ptrdiff_t stop;

	inline ptrdiff_t len() const
	{
		return stop - start;
	}
};

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner MUST guarantee len >= 0
// - Callers MAY choose to handle stride < 0
template<typename T>
struct vec 
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	size_t size() const { return static_cast<size_t>(len); }

	T& operator[](const ptrdiff_t i)
	{
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const
	{
		return ptr[stride * i];
	}

	vec<T> fill(
		const T start = 0,
		const T increment = 0)
	{
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = start + (i * increment);
		return (*this);
	}

	template<typename Tform>
	void transform(Tform op)
	{
		for ( size_t i = 0; i < size(); ++i )
			(*this)[i] = coerce_cast<T>(op((*this)[i]));
	}
	
	template<Unop Op>
	void transform()
	{
		transform(unop<Op,T>{});
	}

	// template<typename Vec, typename Index, typename Reduce>
	// void gather(Vec src, vec<Index> index, Reduce op)
	// {
	// 	assert(size() == index.size());
	// 	for ( size_t i = 0; i < size(); ++i )
	// 		(*this)[i] = 
	//
	vec<T> slice(bounds b) const
	{
		return {
			.ptr = ptr + (stride * b.start), 
			.len = b.len(),
			.stride = stride,
		};
	}
	
	bounds all_elements() const
	{
		return {0, len};
	}
};

// A non-owning strided matrix
// - Owner is responsible for managing memory
// - Owner MUST guarantee nrows >= 0 and ncols >= 0
// - Callers MAY choose to handle row_stride < 0 and col_stride < 0
template<typename T>
struct mat 
{
	T * ptr;
	ptrdiff_t nrows;
	ptrdiff_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;

	vec<T> row(const ptrdiff_t i) const
	{
		return {
			.ptr = ptr + (row_stride * i), 
			.len = ncols, 
			.stride = col_stride
		};
	}

	vec<T> col(const ptrdiff_t i) const
	{
		return {
			.ptr = ptr + (col_stride * i), 
			.len = nrows, 
			.stride = row_stride
		};
	}
	
	mat<T> slice_rows(const bounds b) const
	{
		return {
			.ptr = ptr + (row_stride * b.start),
			.nrows = b.len(),
			.ncols = ncols,
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}

	mat<T> slice_cols(const bounds b) const
	{
		return {
			.ptr = ptr + (col_stride * b.start),
			.nrows = nrows,
			.ncols = b.len(),
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}
	bounds all_rows() const
	{
		return {0, nrows};
	}

	bounds all_cols() const
	{
		return {0, ncols};
	}
};

//// R objects
//--------------
// Initialize struct from R object

#ifdef USING_R
template<typename T>
vec<T> as_vec(SEXP x)
{
	if ( x != R_NilValue )
	{
		return {
			.ptr = data_ptr<T>(x),
			.len = static_cast<ptrdiff_t>(XLENGTH(x)),
			.stride = 1,
		};
	}
	else
	{
		return {nullptr, 0, 0};
	}
}
template<typename T>
mat<T> as_mat(SEXP x)
{
	if ( x != R_NilValue )
	{
		return {
			.ptr = data_ptr<T>(x),
			.nrows = static_cast<ptrdiff_t>(Rf_nrows(x)),
			.ncols = static_cast<ptrdiff_t>(Rf_ncols(x)),
			.row_stride = 1,
			.col_stride = Rf_nrows(x),
		};
	}
	else
	{
		return {nullptr, 0, 0, 0, 0};
	}
}
#endif // USING_R

#endif // CARDINAL_CORE_CORE
