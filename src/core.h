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

// A Num supports arithmetic operations
template<class T>
concept Num = std::is_arithmetic_v<std::remove_cvref_t<T>>;

// A Vec supports 1D tensor operations
// - MUST be trivially copyable as a struct
// - MUST implement ssize() -> ptrdiff_t
// - MUST implement operator[](ptrdiff_t i)
template<class V>
concept Vec = 
	std::is_standard_layout_v<V> &&
	std::is_trivially_copyable_v<V> &&
	requires (const std::remove_cvref_t<V>& v, ptrdiff_t i)
	{
		{ v.ssize() } -> std::convertible_to<ptrdiff_t>;
		{ v[i] } -> Num;
	};

// Proxy type to define UnaryOp and BinaryOp
struct num_arg {
	template<Num T>
	operator T() const;
};

// A callable with one arithmetic argument
template<class F>
concept UnaryOp = std::invocable<F, num_arg>;

// A callable with two arithmetic arguments
template<class F>
concept BinaryOp = std::invocable<F, num_arg, num_arg>;

//// Sentinels
//-------------
// Coercion and NA/missing/incomparable values

// Proxy type to define numeric traits
template<class T>
struct num_traits;

// NaNs are always incomparable so treat them as NAs
template<std::floating_point T>
struct num_traits<T> {
	static constexpr bool is_na(const T x) noexcept {
		return std::isnan(x);
	}
	static constexpr T na_value() noexcept {
		return std::numeric_limits<T>::quiet_NaN();
	}
};

// R NAs
#ifdef USING_R
template<>
struct num_traits<int> {
	static bool is_na(const int x) noexcept { return x == NA_INTEGER; }
	static int na_value() noexcept { return NA_INTEGER; }
};
template<>
struct num_traits<double> {
	static bool is_na(const double x) noexcept { return std::isnan(x); }
	static double na_value() noexcept { return NA_REAL; }
};
#endif // USING_R

// A HasNA type supports NA sentinel values
template<class T>
concept HasNA = Num<T> && requires (const T& x)
	{
		{ num_traits<T>::is_na(x) } -> std::same_as<bool>;
		{ num_traits<T>::na_value() } -> std::same_as<T>;
	};

// NA values (lowest if not defined so na_value<bool> -> false)
template<Num T>
constexpr T na_value() noexcept
{
	if constexpr ( HasNA<T> )
		return num_traits<T>::na_value();
	else
		return std::numeric_limits<T>::lowest();
}

// Check if a value is NA/missing/incomparable
template<Num T>
constexpr bool is_na(const T x) noexcept
{
	if constexpr ( HasNA<T> )
		return num_traits<T>::is_na(x);
	else
		return false;
}

// Coerce while preserving NAs across types if possible
template<Num Out, Num In>
constexpr Out coerce_cast(In x) noexcept
{
	if constexpr ( std::is_same_v<Out,In> )
		return x;
	else
	{
		if ( is_na(x) )
			return na_value<Out>();
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

//// Indexing
//-----------
// Indexes and bounds

// Half-open [start, stop) index bounds
struct bounds
{
	ptrdiff_t start;
	ptrdiff_t stop;

	inline ptrdiff_t width() const noexcept
	{
		return stop - start;
	}
};

// Multidimensional index
template<int N>
struct loc
{
	ptrdiff_t index[N];
};

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
	if ( is_na(x) )
		return na_value<T>();
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
	if ( is_na(lhs) || is_na(rhs) )
		return na_value<T>();
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
			return na_value<Out>();
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
		if ( is_na(ii) )
			return na_value<T>();
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

//// Vector operators
//---------------------
// Deferred math and arithmetic

// Universal unary functions for Vecs
template<Unop Op, Num T = double, Vec V>
Vec auto ufunc(V x) noexcept
{
	return transform<Op,T>(x);
}

// Universal binary functions for Vecs
template<Binop Op, Num T = double, Vec L, Vec R>
Vec auto ufunc(L lhs, R rhs) noexcept
{
	return transform<Op,T>(lhs, rhs);
}

// Vec + Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator+(L lhs, R rhs) noexcept
{
	return ufunc<Add,T>(lhs, rhs);
}

// Vec - Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator-(L lhs, R rhs) noexcept
{
	return ufunc<Subtract,T>(lhs, rhs);
}

// Vec * Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator*(L lhs, R rhs) noexcept
{
	return ufunc<Multiply,T>(lhs, rhs);
}

// Vec / Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator/(L lhs, R rhs) noexcept
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
	T value;
	ptrdiff_t len;

	constexpr ptrdiff_t ssize() const noexcept
	{
		return len;
	}

	constexpr T operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return value;
	}
};

// Generator vector yielding a sequence
template<Num T = double>
struct seq
{
	T start;
	ptrdiff_t len;
	T step = 1;

	constexpr ptrdiff_t ssize() const noexcept
	{
		return len;
	}

	constexpr T operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return start + (i * step);
	}
};

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner is responsible for data validity
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

	vec<T>& fill(const T value = 0) noexcept
	{
		return this->assign(rep<T>{value, len});
	}

	vec<T>& seqfill(const T start = 0, const T step = 1) noexcept
	{
		return this->assign(seq<T>{start, len, step});
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
	template<Binop Op, Vec V, BinaryOp Tform = binop<Op,T>>
	vec<T>& transform(const V src, const Tform op = Tform{}) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		return (*this);
	}
	// Assign (*this)[i] = src[index[i]] for i in index
	template<Binop Op = Rhs, Vec Index, Vec V>
	vec<T>& gather(const Index index, const V src) noexcept
	{
		assert(index.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			auto ii = index[i];
			if ( is_na(ii) )
			{
				(*this)[i] = na_value<T>;
				continue;
			}
			(*this)[i] = ufunc<Op,T>((*this)[i], coerce_cast<T>(src[ii]));
		}
		return (*this);
	}

	// Assign (*this)[index[i]] = src[i] for i in index
	template<Binop Op = Rhs, Vec Index, Vec V>
	vec<T>& scatter(const Index index, const V src) noexcept
	{
		assert(index.ssize() == src.ssize());
		for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
		{
			auto ii = index[i];
			if ( is_na(ii) )
				continue;
			(*this)[ii] = ufunc<Op,T>((*this)[ii], coerce_cast<T>(src[i]));
		}
		return (*this);
	}

	// vec<T> += Vec
	template<Vec V>
	vec<T>& operator+=(const V src) noexcept
	{
		return this->transform<Add>(src);
	}
	
	// vec<T> -= Vec
	template<Vec V>
	vec<T>& operator-=(const V src) noexcept
	{
		return this->transform<Subtract>(src);
	}

	// vec<T> *= Vec
	template<Vec V>
	vec<T>& operator*=(const V src) noexcept
	{
		return this->transform<Multiply>(src);
	}

	// vec<T> /= Vec
	template<Vec V>
	vec<T>& operator/=(const V src) noexcept
	{
		return this->transform<Divide>(src);
	}

	// Return a sliced view from b.start to b.stop
	vec<T> slice(const bounds b) const noexcept
	{
		return {
			.ptr = ptr + (stride * b.start), 
			.len = b.width(),
			.stride = stride,
		};
	}
	
	bounds all_elements() const noexcept
	{
		return {0, len};
	}

	#ifdef USING_R
	static vec<T> from(SEXP x) noexcept
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
	#endif // USING_R
};

// An locally-scoped owning vector
template<Num T>
struct local_vec : vec<T>
{
	explicit local_vec(ptrdiff_t n) :
		vec<T>{
			.ptr = new T[n],
			.len = n,
			.stride = 1,
		} {}
	~local_vec() { delete[] this->ptr; }

	local_vec(const local_vec&) = delete;
	local_vec(local_vec&&) = delete;
	local_vec& operator=(const local_vec&) = delete;
	local_vec& operator=(local_vec&&) = delete;

	vec<T> borrow() const noexcept
	{
		return {
			.ptr = this->ptr,
			.len = this->len,
			.stride = 1,
		};
	}
};

//// Matrices
//------------
// 2D array operations

// A non-owning strided matrix
// - Owner is responsible for managing memory
// - Owner is responsible for data validity
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
			.nrows = b.width(),
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
			.ncols = b.width(),
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

	#ifdef USING_R
	static mat<T> from(SEXP x) noexcept
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
};

#endif // CARDINAL_CORE_CORE
