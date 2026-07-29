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
constexpr bool is_incomparable(T x)
{
	if constexpr ( std::is_floating_point_v<T> )
		return std::isnan(x);
	else
		return x == make_incomparable<T>();
}

#ifdef USING_R
template<>
inline bool is_incomparable<int>(int x) { return x == NA_INTEGER; }
template<>
inline bool is_incomparable<double>(double x) { return ISNAN(x); }
#endif // USING_R

template<typename Out, typename In>
Out coerce_cast(In x)
{
	if ( is_incomparable(x) )
		return make_incomparable<Out>();
	else
		return static_cast<Out>(x);
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
// Universal unary functions

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
T ufunc(T x)
{
	if constexpr ( Op != Identity )
	{
		if ( is_incomparable(x) )
			return make_incomparable<T>();
	}
	if constexpr ( Op == Identity )
		return x;
	else if constexpr ( Op == Abs )
		return std::abs(x);
	else if constexpr ( Op == Log )
		return std::log(x);
	else if constexpr ( Op == Log2 )
		return std::log2(x);
	else if constexpr ( Op == Log1p )
		return std::log1p(x);
	else if constexpr ( Op == Exp )
		return std::exp(x);
	else if constexpr ( Op == Exp2)
		return std::exp2(x);
	else if constexpr ( Op == Expm1 )
		return std::expm1(x);
	else
		static_assert(Op != Op, "ufunc: unsupported unary op");
}

template<Unop Op, typename Out = double, typename In = Out>
struct unop {
	Out operator()(In x) const
	{
		if constexpr ( std::is_same_v<Out,In> )
			return ufunc<Op,Out>(x);
		else
			return coerce_cast<Out>(ufunc<Op,In>(x));
	}
};

//// Binary operations
//---------------------
// Universal binary functions

enum Binop {
	Lhs,
	Rhs,
	Add,
	Subtract,
	Multiply,
	Divide,
	Max,
	Min
};

template<Binop Op, typename T = double>
T ufunc(T lhs, T rhs)
{
	if constexpr ( Op != Lhs && Op != Rhs )
	{
		if ( is_incomparable(lhs) || is_incomparable(rhs) )
			return make_incomparable<T>();
	}
	if constexpr ( Op == Lhs )
		return lhs;
	else if constexpr ( Op == Rhs )
		return rhs;
	else if constexpr ( Op == Add )
		return lhs + rhs;
	else if constexpr ( Op == Subtract )
		return lhs - rhs;
	else if constexpr ( Op == Multiply )
		return lhs * rhs;
	else if constexpr ( Op == Divide )
		return lhs / rhs;
	else if constexpr ( Op == Max )
		return lhs > rhs ? lhs : rhs;
	else if constexpr ( Op == Min )
		return lhs < rhs ? lhs : rhs;
	else
		static_assert(Op != Op, "ufunc: unsupported binary op");
}

template<Binop Op, typename Out = double, typename In = Out>
struct binop {
	static Out identity()
	{
		if constexpr ( Op == Add )
			return 0;
		else if constexpr ( Op == Subtract )
			return 0;
		else if constexpr ( Op == Multiply )
			return 1;
		else if constexpr ( Op == Divide )
			return 1;
		else if constexpr ( Op == Max )
		{
			if constexpr ( std::numeric_limits<Out>::has_infinity )
				return -std::numeric_limits<Out>::infinity();
			else
				return std::numeric_limits<Out>::lowest();
		}
		else if constexpr ( Op == Min )
		{
			if constexpr ( std::numeric_limits<Out>::has_infinity )
				return std::numeric_limits<Out>::infinity();
			else
				return std::numeric_limits<Out>::max();
		}
		else
			return make_incomparable<Out>();
	}
	Out operator()(In lhs, In rhs) const
	{
		if constexpr ( std::is_same_v<Out,In> )
			return ufunc<Op,Out>(lhs, rhs);
		else
			return coerce_cast<Out>(ufunc<Op,In>(lhs, rhs));
	}
};

//// Vectors
//-----------
// 1D array operations

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
template<typename T>
struct vec 
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	size_t size() const 
	{
		return static_cast<size_t>(len);
	}

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

	template<Unop Op, typename Tform = unop<Op,T>>
	vec<T> transform(Tform op = Tform{})
	{
		for ( size_t i = 0; i < size(); ++i )
			(*this)[i] = op((*this)[i]);
		return (*this);
	}
	
	template<Binop Op, typename Vec, typename Tform = binop<Op,T>>
	vec<T> transform(Vec src, Tform op = Tform{})
	{
		assert(this->size() == src.size());
		for ( size_t i = 0; i < size(); ++i )
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}

	template<Binop Op = Rhs, typename Index, typename Vec>
	vec<T> gather(vec<Index> index, Vec src)
	{
		assert(this->size() == index.size());
		binop<Op,T> op{};
		for ( size_t i = 0; i < size(); ++i )
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[index[i]]));
		return (*this);
	}

	template<Binop Op = Rhs, typename Index, typename Vec>
	vec<T> scatter(vec<Index> index, Vec src)
	{
		assert(src.size() == index.size());
		binop<Op,T> op{};
		for ( size_t i = 0; i < index.size(); ++i )
			(*this)[index[i]] = op((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}

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

template<typename Vec, typename Reduce, typename T = double>
T reduce(Vec input, Reduce op, T init) 
{
	T output = init;
	for ( size_t i = 0; i < input.size(); ++i )
		output = op(output, coerce_cast<T>(input[i]));
	return output;
}

template<Binop Op, typename Vec, typename T = double>
T reduce(Vec input) 
{
	return reduce(input, binop<Op,T>{}, binop<Op,T>::identity());
}

//// Matrices
//------------
// 2D array operations

// A non-owning strided matrix
// - Owner is responsible for managing memory
// - Owner MUST guarantee nrows >= 0 and ncols >= 0
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
