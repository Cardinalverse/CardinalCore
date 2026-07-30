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

template<class T>
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

template<class T>
constexpr bool is_incomparable(T x)
{
	if constexpr ( std::is_floating_point_v<T> )
		return std::isnan(x);
	else
		return x == make_incomparable<T>();
}

template<class Out, class In>
Out coerce_cast(In x)
{
	if constexpr ( std::is_same_v<Out,In> )
		return x;
	else
	{
		if ( is_incomparable(x) )
			return make_incomparable<Out>();
		else
			return static_cast<Out>(x);
	}
}

//// Data Pointer
//----------------
// Get a mutable data pointer from a managed object

#ifdef USING_R
template<class T>
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

template<Unop Op, class T = double>
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

template<Unop Op, class Out = double, class In = Out>
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

template<Binop Op, class T = double>
T ufunc(T lhs, T rhs)
{
	if constexpr ( Op != Lhs && Op != Rhs )
	{
		if ( is_incomparable(lhs) || is_incomparable(rhs) )
			return make_incomparable<T>();
	}
	if constexpr ( Op == Lhs || Op == Rhs )
	{
		if constexpr ( Op == Lhs )
			return lhs;
		if constexpr ( Op == Rhs )
			return rhs;
	}
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

template<Binop Op, class Out = double, class In = Out>
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

//// Vector expressions
//---------------------
// Generic operations on vectors

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

// Vector with elementwise unary transformation
// - MUST implement (ptrdiff_t)ssize() and (T)operator()(ptrdiff_t)
template<class Vec, class Tform, class T = double>
struct vec_unop
{
	Vec x;
	Tform op;

	ptrdiff_t ssize() const
	{
		return x.ssize();
	}

	T operator[](ptrdiff_t i) const
	{
		return coerce_cast<T>(op(x[i]));
	}
};

// Vector with elementwise binary transformation
// - MUST implement (ptrdiff_t) ssize() and (T) operator()(ptrdiff_t)
template<class LVec, class RVec, class Tform, class T = double>
struct vec_binop
{
	LVec lhs;
	RVec rhs;
	Tform op;

	ptrdiff_t ssize() const
	{
		return lhs.ssize();
	}

	T operator[](ptrdiff_t i) const
	{
		return coerce_cast<T>(op(lhs[i], rhs[i]));
	}
};

// Vector subscripted at the given indices
// - MUST implement (ptrdiff_t)ssize() and (T)operator()(ptrdiff_t)
template<class Vec, class Index, class T = double>
struct vec_indexed
{
	Vec x;
	Index index;

	ptrdiff_t ssize() const
	{
		return index.ssize();
	}

	T operator[](ptrdiff_t i) const
	{
		return coerce_cast<T>(x[index[i]]);
	}
};

// Check for incomparables
template<class Vec>
bool any_incomparable(const Vec input)
{
	for ( ptrdiff_t i = 0; i < input.ssize(); ++i )
		if ( is_incomparable(input[i]) )
			return true;
	return false;
}

// Transform with elementwise unary functor
template<Unop Op, 
	class T = double, 
	class Vec, 
	class Tform = unop<Op,T>
	>
vec_unop<Vec,Tform,T> transform(
	const Vec input, 
	const Tform op = Tform{})
{
	return {input, op};
}

// Transform with elementwise binary functor
template<
	Binop Op, 
	class T = double, 
	class LVec, 
	class RVec, 
	class Tform = binop<Op,T>
	>
vec_binop<LVec,RVec,Tform,T> transform(
	const LVec lhs, 
	const RVec rhs,
	const Tform op = Tform{})
{
	assert(lhs.ssize() == rhs.ssize());
	return {lhs, rhs, op};
}

// Gather vector elements and given indices
template<
	class T = double,
	class Index,
	class Vec
	>
vec_indexed<Vec,Index,T> gather(
	const Index index,
	const Vec input)
{
	return {input, index};
}

// Reduce vector elements with binary functor
template<
	Binop Op, 
	class T = double, 
	class Vec, 
	class Reduce = binop<Op,T>
	>
T reduce(
	const Vec input,
	const Reduce op = Reduce{},
	const T init = Reduce::identity()) 
{
	T output = init;
	for ( ptrdiff_t i = 0; i < input.ssize(); ++i )
		output = op(output, coerce_cast<T>(input[i]));
	return output;
}

//// Vectors
//-----------
// 1D array operations

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner MUST guarantee len >= 0
template<class T>
struct vec 
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	ptrdiff_t ssize() const
	{
		return len;
	}

	T& operator[](const ptrdiff_t i)
	{
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const
	{
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	vec<T> fill(const T start = 0, const T increment = 0)
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			(*this)[i] = start + (i * increment);
		return (*this);
	}

	// Elementwise assignment
	template<class Vec>
	vec<T> assign(const Vec src)
	{
		assert(src.ssize() == this->ssize());
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			(*this)[i] = coerce_cast<T>(src[i]);
		return (*this);
	}

	// Elementwise in-place unary transformations
	template<Unop Op, class Tform = unop<Op,T>>
	vec<T> transform(const Tform op = Tform{})
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			(*this)[i] = op((*this)[i]);
		return (*this);
	}
	
	// Elementwise in-place binary transformations
	template<Binop Op, class Vec, class Tform = binop<Op,T>>
	vec<T> transform(const Vec src, const Tform op = Tform{})
	{
		assert(src.ssize() == this->ssize());
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}

	// Assign (*this)[i] = src[index[i]] for i in index
	template<Binop Op = Rhs, class Index, class Vec>
	vec<T> gather(const vec<Index> index, const Vec src)
	{
		assert(index.ssize() == this->ssize());
		for ( ptrdiff_t i = 0; i < index.ssize(); ++i )
			(*this)[i] = ufunc<Op,T>((*this)[i], coerce_cast<T>(src[index[i]]));
		return (*this);
	}

	// Assign (*this)[index[i]] = src[i] for i in index
	template<Binop Op = Rhs, class Index, class Vec>
	vec<T> scatter(const vec<Index> index, const Vec src)
	{
		assert(index.ssize() == src.ssize());
		for ( ptrdiff_t i = 0; i < index.ssize(); ++i )
			(*this)[index[i]] = ufunc<Op,T>((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}

	// Return a sliced view from b.start to b.stop
	vec<T> slice(const bounds b) const
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

//// Matrices
//------------
// 2D array operations

// A non-owning strided matrix
// - Owner is responsible for managing memory
// - Owner MUST guarantee nrows >= 0 and ncols >= 0
template<class T>
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
template<class T>
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
template<class T>
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
