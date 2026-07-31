#ifndef CARDINAL_CORE_CORE
#define CARDINAL_CORE_CORE

#include <cassert>
#include <cstddef>
#include <cmath>
#include <limits>
#include <concepts>
#include <type_traits>

using size_t = std::size_t;
using ptrdiff_t = std::ptrdiff_t;

//// Concepts
//------------
// Traits supporting template dispatch

// A Vec supports 1D tensor operations
// - MUST be trivially copyable as a struct
// - MUST implement ssize() -> ptrdiff_t
// - MUST implement operator[](ptrdiff_t i)
template<class V>
concept Vec = 
	std::is_trivial_v<V> &&
	std::is_standard_layout_v<V> &&
	requires (const std::remove_cvref_t<V>& v, ptrdiff_t i) {
		{ v.ssize() } -> std::convertible_to<ptrdiff_t>;
		{ v[i] };
	};

// A Num supports arithmetic operations
template<class T>
concept Num = std::is_arithmetic_v<T>;

// Proxy type to define UnaryOp and BinaryOp
struct num_arg {
	template<Num T>
	operator T() const;
};

// A callable with one arithmetic argument
template<class Tform>
concept UnaryOp = std::invocable<Tform, num_arg>;

// A callable with two arithmetic arguments
template<class Tform>
concept BinaryOp = std::invocable<Tform, num_arg, num_arg>;

//// Data traits
//---------------
// Coercion and incomparable (missing) values

// Proxy type to define incomparables
template<Num T>
struct num_traits;

// NaNs are always incomparable
template<Num T>
struct num_traits {
	static constexpr T incomparable() noexcept {
		return std::numeric_limits<T>::quiet_NaN();
	}
};

// R NAs
#ifdef USING_R
template<>
struct num_traits<int> {
	static int incomparable() noexcept { return NA_INTEGER; }
};
template<>
struct num_traits<double> {
	static double incomparable() noexcept { return NA_REAL; }
};
#endif // USING_R

// Incomparable values
template<Num T>
constexpr T incomparable() noexcept {
	return num_traits<T>::incomparable();
}

// A MaybeIncomparble type might be incomparable
template<class T>
concept MaybeIncomparable = requires { num_traits<T>::incomparable(); };

// Check if a value is comparable
template<Num T>
constexpr bool is_incomparable(const T x) noexcept
{
	if constexpr ( std::is_floating_point_v<T> )
		return std::isnan(x);
	else if constexpr ( MaybeIncomparable<T> )
		return x == incomparable<T>();
	else
		return false;
}

// Coerce while preserving incomparables across types
template<Num Out, Num In>
constexpr Out coerce_cast(In x) noexcept
{
	if constexpr ( std::is_same_v<Out,In> )
		return x;
	else
	{
		if ( is_incomparable(x) )
			return incomparable<Out>();
		else
			return static_cast<Out>(x);
	}
}

//// Data Pointer
//----------------
// Get a mutable data pointer from a managed object

#ifdef USING_R
template<class T>
T * data_ptr(SEXP x) noexcept;
template<>
inline int * data_ptr<int>(SEXP x) noexcept { return INTEGER(x); }
template<>
inline double * data_ptr<double>(SEXP x) noexcept { return REAL(x); }
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

template<Unop Op, Num T = double>
constexpr T ufunc(T x) noexcept
{

	if constexpr ( Op == Identity )
		return x;
	if ( is_incomparable(x) )
		return incomparable<T>();
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

template<Unop Op, Num Out = double, Num In = Out>
struct unop {
	Out operator()(In x) const noexcept
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

template<Binop Op, Num T = double>
constexpr T ufunc(T lhs, T rhs) noexcept
{
	if constexpr ( Op == Lhs )
		return lhs;
	if constexpr ( Op == Rhs )
		return rhs;
	if ( is_incomparable(lhs) || is_incomparable(rhs) )
		return incomparable<T>();
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

template<Binop Op, Num Out = double, Num In = Out>
struct binop {
	static Out identity() noexcept
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
				return +std::numeric_limits<Out>::infinity();
			else
				return std::numeric_limits<Out>::max();
		}
		else
			return incomparable<Out>();
	}
	Out operator()(In lhs, In rhs) const noexcept
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

// Half-open [start, stop) index bounds
struct bounds 
{
	ptrdiff_t start;
	ptrdiff_t stop;

	inline ptrdiff_t len() const noexcept
	{
		return stop - start;
	}
};

// Vector subscripted at the given indices
template<Vec V, Vec Index, Num T = double>
struct vec_indexed
{
	V x;
	Index index;

	ptrdiff_t ssize() const noexcept
	{
		return index.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		auto ii = index[i];
		if ( is_incomparable(ii) )
			return incomparable<T>();
		else
			return coerce_cast<T>(x[ii]);
	}
};

// Gather vector elements and given indices
template<Num T = double, Vec Index, Vec V>
constexpr auto gather(
	const Index index,
	const V input) noexcept -> vec_indexed<V,Index,T>
{
	return {
		.x = input, 
		.index = index,
	};
}

// Vector with elementwise unary transformation
template<Vec V, UnaryOp Tform, Num T = double>
struct vec_unop
{
	V x;
	Tform op;

	ptrdiff_t ssize() const noexcept
	{
		return x.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		return coerce_cast<T>(op(x[i]));
	}
};

// Transform with elementwise unary functor
template<Unop Op, Num T = double, Vec V, UnaryOp Tform = unop<Op,T>>
constexpr auto transform(
	const V input, 
	const Tform op = Tform{}) noexcept -> vec_unop<V,Tform,T>
{
	return {
		.x = input, 
		.op = op,
	};
}

// Vector with elementwise binary transformation
template<Vec L, Vec R, BinaryOp Tform, Num T = double>
struct vec_binop
{
	L lhs;
	R rhs;
	Tform op;

	ptrdiff_t ssize() const noexcept
	{
		return lhs.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		return coerce_cast<T>(op(lhs[i], rhs[i]));
	}
};

// Transform with elementwise binary functor
template<Binop Op, Num T = double, Vec L, Vec R, BinaryOp Tform = binop<Op,T>>
constexpr auto transform(
	const L lhs, 
	const R rhs,
	const Tform op = Tform{}) noexcept -> vec_binop<L,R,Tform,T>
{
	assert(lhs.ssize() == rhs.ssize());
	return {
		.lhs = lhs, 
		.rhs = rhs, 
		.op = op,
	};
}

// Reduce vector elements with binary functor
template<Binop Op, Num T = double, Vec V, BinaryOp Reduce = binop<Op,T>>
T reduce(
	const V input,
	const Reduce op = Reduce{},
	const T init = Reduce::identity()) noexcept
{
	T output = init;
	for ( ptrdiff_t i = 0; i < input.ssize(); ++i )
		output = op(output, coerce_cast<T>(input[i]));
	return output;
}

// Universal unary functions for Vecs
template<Unop Op, Num T = double, Vec V>
auto ufunc(V x)
{
	return transform<Op,T>(x);
}

// Universal binary functions for Vecs
template<Binop Op, Num T = double, Vec L, Vec R>
auto ufunc(L lhs, R rhs)
{
	return transform<Op,T>(lhs, rhs);
}

//// Vector operators
//---------------------
// Deferred vectorized arithmetic

// Vec + Vec
template<Num T = double, Vec L, Vec R>
constexpr auto operator+(L lhs, R rhs) noexcept
{
	return ufunc<Add,T>(lhs, rhs);
}

// Vec - Vec
template<Num T = double, Vec L, Vec R>
constexpr auto operator-(L lhs, R rhs) noexcept
{
	return ufunc<Subtract,T>(lhs, rhs);
}

// Vec * Vec
template<Num T = double, Vec L, Vec R>
constexpr auto operator*(L lhs, R rhs) noexcept
{
	return ufunc<Multiply,T>(lhs, rhs);
}

// Vec / Vec
template<Num T = double, Vec L, Vec R>
constexpr auto operator/(L lhs, R rhs) noexcept
{
	return ufunc<Divide,T>(lhs, rhs);
}

//// Vectors
//-----------
// 1D array operations

// Generator vector repeating a constant
template<Num T = double>
struct rep
{
	T x;
	ptrdiff_t len;

	constexpr ptrdiff_t ssize() const noexcept
	{
		return len;
	}

	constexpr T operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return x;
	}
};

// Generator vector yielding a sequence
template<Num T = double>
struct seq
{
	T start;
	ptrdiff_t len;
	T increment = 1;

	constexpr ptrdiff_t ssize() const noexcept
	{
		return len;
	}

	constexpr T operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return start + (i * increment);
	}
};

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner MUST guarantee len >= 0
template<Num T>
struct vec
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	ptrdiff_t ssize() const noexcept
	{
		return len;
	}

	T& operator[](const ptrdiff_t i) noexcept
	{
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	vec<T>& fill(const T start = 0, const T increment = 0) noexcept
	{
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = start + (i * increment);
		return (*this);
	}

	vec<T>& seq_fill() noexcept
	{
		return this->fill(0, 1);
	}

	// Elementwise assignment
	template<Vec V>
	vec<T>& assign(const V src) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = coerce_cast<T>(src[i]);
		return (*this);
	}

	// Elementwise in-place unary transformations
	template<Unop Op, UnaryOp Tform = unop<Op,T>>
	vec<T>& transform(
		const Tform op = Tform{}) noexcept
	{
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = op((*this)[i]);
		return (*this);
	}

	// Elementwise in-place binary transformations
	template<Binop Op, Vec Src, BinaryOp Tform = binop<Op,T>>
	vec<T>& transform(const Src src, const Tform op = Tform{}) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}
	// Assign (*this)[i] = src[index[i]] for i in index
	template<Binop Op = Rhs, Vec Index, Vec Src>
	vec<T>& gather(const Index index, const Src src) noexcept
	{
		assert(index.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			auto ii = index[i];
			if ( is_incomparable(ii) )
			{
				(*this)[i] = incomparable<T>;
				continue;
			}
			(*this)[i] = ufunc<Op,T>((*this)[i], coerce_cast<T>(src[ii]));
		}
		return (*this);
	}

	// Assign (*this)[index[i]] = src[i] for i in index
	template<Binop Op = Rhs, Vec Index, Vec Src>
	vec<T>& scatter(const Index index, const Src src) noexcept
	{
		assert(index.ssize() == src.ssize());
		for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
		{
			auto ii = index[i];
			if ( is_incomparable(ii) )
				continue;
			(*this)[ii] = ufunc<Op,T>((*this)[ii], coerce_cast<T>(src[i]));
		}
		return (*this);
	}

	// vec + vec
	template<Vec Src>
	vec<T>& operator+=(const Src src) noexcept
	{
		return this->transform<Add>(src);
	}
	
	// vec - vec
	template<Vec Src>
	vec<T>& operator-=(const Src src) noexcept
	{
		return this->transform<Subtract>(src);
	}

	// vec * vec
	template<Vec Src>
	vec<T>& operator*=(const Src src) noexcept
	{
		return this->transform<Multiply>(src);
	}

	// vec / vec
	template<Vec Src>
	vec<T>& operator/=(const Src src) noexcept
	{
		return this->transform<Divide>(src);
	}

	// Return a sliced view from b.start to b.stop
	vec<T> slice(const bounds b) const noexcept
	{
		return {
			.ptr = ptr + (stride * b.start), 
			.len = b.len(),
			.stride = stride,
		};
	}
	
	bounds all_elements() const noexcept
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
template<Num T>
struct mat
{
	T * ptr;
	ptrdiff_t nrows;
	ptrdiff_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;

	vec<T> row(const ptrdiff_t i) const noexcept
	{
		return {
			.ptr = ptr + (row_stride * i), 
			.len = ncols, 
			.stride = col_stride,
		};
	}

	vec<T> col(const ptrdiff_t i) const noexcept
	{
		return {
			.ptr = ptr + (col_stride * i), 
			.len = nrows, 
			.stride = row_stride,
		};
	}
	
	mat<T> slice_rows(const bounds b) const noexcept
	{
		return {
			.ptr = ptr + (row_stride * b.start),
			.nrows = b.len(),
			.ncols = ncols,
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}

	mat<T> slice_cols(const bounds b) const noexcept
	{
		return {
			.ptr = ptr + (col_stride * b.start),
			.nrows = nrows,
			.ncols = b.len(),
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}
	bounds all_rows() const noexcept
	{
		return {0, nrows};
	}

	bounds all_cols() const noexcept
	{
		return {0, ncols};
	}
};

//// R objects
//--------------
// Initialize struct from R object

#ifdef USING_R
template<class T>
vec<T> as_vec(SEXP x) noexcept
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
mat<T> as_mat(SEXP x) noexcept
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
