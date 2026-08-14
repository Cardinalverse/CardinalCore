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

// A Masked Vec carries a mask for valid/invalid elements
template<class M>
concept Masked = Vec<M> &&
	requires(const std::remove_cvref_t<M>& m, ptrdiff_t i)
	{
		{ m.is_valid(i) } -> std::same_as<bool>;
		{ m.get_data() } -> Vec;
		{ m.get_mask() } -> Vec;
	};

// Check if a Vec element is valid
template<Vec V>
constexpr bool is_valid(V& x, ptrdiff_t i)
{
	if constexpr ( Masked<V> )
		return x.is_valid(i);
	else
		return true;
}

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

// A HasNA type supports NA sentinel values
template<class T>
concept HasNA = Num<T> && requires (const T& x)
	{
		{ num_traits<T>::na_value() } -> std::same_as<T>;
		{ num_traits<T>::is_na(x) } -> std::same_as<bool>;
	};

// NA values (defaults to numeric_limits<T>::lowest() if undefined)
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

// R NAs
#ifdef USING_R
template<>
struct num_traits<int> {
	static int na_value() noexcept { return NA_INTEGER; }
	static bool is_na(const int x) noexcept { return x == NA_INTEGER; }
};
template<>
struct num_traits<double> {
	static double na_value() noexcept { return NA_REAL; }
	static bool is_na(const double x) noexcept { return std::isnan(x); }
};
#endif // USING_R

// Count of invalid/missing/NA items in x
template<Vec V>
ptrdiff_t n_missing(V x) noexcept
{
	ptrdiff_t count = 0;
	for ( ptrdiff_t i = 0; i < x.ssize(); ++i )
	{
		if constexpr ( Masked<V> )
			count += !x.is_valid(i) || is_na(x[i]);
		else
			count += is_na(x[i]);
	}
	return count;
}

// Count of valid/non-missing/non-NA items in x
template<Vec V>
ptrdiff_t n_present(V x) noexcept
{
	return x.ssize() - n_missing(x);
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
};

// Multidimensional index
template<int N>
struct loc
{
	ptrdiff_t loc[N];

	ptrdiff_t operator[](ptrdiff_t i) const noexcept
	{
		return loc[i];
	}
};

//// Unary operations
//--------------------
// Universal unary functions

enum Unop {
	// Identity
	Identity, IsNA, NotNA,
	// Logic
	Not,
	// Math
	Abs, Sign, Log, Log2, Log1p, Exp, Exp2, Expm1,
};

template<Unop Op, Num T = double>
constexpr T ufunc(T x) noexcept
{
	// Identity
	if constexpr ( Op == Identity )
		return x;
	else if constexpr ( Op == IsNA )
		return is_na(x);
	else if constexpr ( Op == NotNA )
		return !is_na(x);
	// Logic
	else if constexpr ( Op == Not )
		return is_na(x) ? coerce_cast<T>(x) : !x;
	else
	{
		// NAs
		if ( is_na(x) )
			return na_value<T>();
		// Math
		else if constexpr ( Op == Abs )
			return std::abs(x);
		else if constexpr ( Op == Sign )
			return x ? std::copysign(1, x) : 0;
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
		// Not implemented
		else
			static_assert(Op != Op, "ufunc: unsupported unary op");
	}
}

template<Unop Op, Num T = double, Vec V>
constexpr Vec auto ufunc(V x) noexcept;

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
	// Identity
	Lhs, Rhs,
	// Logic
	And, Or,
	// Compare
	Eq, Ne, Lt, Le, Gt, Ge,
	// Arithmetic
	Add, Sub, Mul, Div, Max, Min
};

template<Binop Op, Num T = double>
constexpr T ufunc(T lhs, T rhs) noexcept
{
	// Identity
	if constexpr ( Op == Lhs )
		return lhs;
	else if constexpr ( Op == Rhs )
		return rhs;
	// Logic
	else if constexpr ( Op == And )
	{
		if ( is_na(lhs) && is_na(rhs) )
			return na_value<T>();
		else if ( is_na(lhs) )
			return rhs ? na_value<T>() : false;
		else if ( is_na(rhs) )
			return lhs ? na_value<T>() : false;
		else
			return lhs && rhs;
	}
	else if constexpr ( Op == Or )
	{
		if ( is_na(lhs) && is_na(rhs) )
			return na_value<T>();
		else if ( is_na(lhs) )
			return rhs ? true : na_value<T>();
		else if ( is_na(rhs) )
			return lhs ? true : na_value<T>();
		else
			return lhs || rhs;
	}
	else
	{
		// NAs
		// if ( is_na(lhs) || is_na(rhs) )
		// 	return na_value<T>();
		// Compare
		// else if constexpr ( Op == Eq )
		if constexpr ( Op == Eq )
			return lhs == rhs;
		else if constexpr ( Op == Ne )
			return lhs != rhs;
		else if constexpr ( Op == Lt )
			return lhs < rhs;
		else if constexpr ( Op == Le )
			return lhs <= rhs;
		else if constexpr ( Op == Gt )
			return lhs > rhs;
		else if constexpr ( Op == Ge )
			return lhs >= rhs;
		// Arithmetic
		else if constexpr ( Op == Add )
			return lhs + rhs;
		else if constexpr ( Op == Sub )
			return lhs - rhs;
		else if constexpr ( Op == Mul )
			return lhs * rhs;
		else if constexpr ( Op == Div )
			return lhs / rhs;
		else if constexpr ( Op == Max )
			return lhs > rhs ? lhs : rhs;
		else if constexpr ( Op == Min )
			return lhs < rhs ? lhs : rhs;
		// Not implemented
		else
			static_assert(Op != Op, "ufunc: unsupported binary op");
	}
}

template<Binop Op, Num T = double, Vec L, Vec R>
constexpr Vec auto ufunc(L lhs, R rhs) noexcept;

template<Binop Op, Num Out = double, Num In = Out>
struct binop {
	static Out identity() noexcept
	{
		if constexpr ( Op == Add )
			return 0;
		else if constexpr ( Op == Sub )
			return 0;
		else if constexpr ( Op == Mul )
			return 1;
		else if constexpr ( Op == Div )
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

// Vector with validity mask
template<Vec V, Vec Mask, Num T = double>
struct vec_masked
{
	V data;
	Mask mask;

	ptrdiff_t ssize() const noexcept
	{
		return data.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		return is_valid(i) ? coerce_cast<T>(data[i]) : na_value<T>();
	}

	bool is_valid(ptrdiff_t i) const noexcept
	{
		return coerce_cast<bool>(mask[i]);
	}

	constexpr V get_data() const noexcept
	{
		return data;
	}

	constexpr Mask get_mask() const noexcept
	{
		return mask;
	}
};

// Mask join
template<Vec L, Vec R, Num T = bool>
struct vec_valid
{
	L ml;
	R mr;

	ptrdiff_t ssize() const noexcept
	{
		return ml.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		if ( is_valid(ml, i) && is_valid(mr, i) )
			return coerce_cast<bool>(ml) && coerce_cast<bool>(mr);
		else
			return false;
	}
};

// Mask with ternary AND logic
template<Vec L, Vec R, Num T = bool>
struct vec_valid_and
{
	L lhs;
	R rhs;

	ptrdiff_t ssize() const noexcept
	{
		return lhs.ssize();
	}

	// Valid if both valid or if either is valid-FALSE
	T operator[](ptrdiff_t i) const noexcept
	{
		bool vl = is_valid(lhs, i);
		bool vr = is_valid(rhs, i);
		return (vl && vr) || (vl && !lhs[i]) || (vr && !rhs[i]);
	}
};

// Mask with ternary OR logic
template<Vec L, Vec R, Num T = bool>
struct vec_valid_or
{
	L lhs;
	R rhs;

	ptrdiff_t ssize() const noexcept
	{
		return lhs.ssize();
	}

	// Valid if both valid or if either is valid-TRUE
	T operator[](ptrdiff_t i) const noexcept
	{
		bool vl = is_valid(lhs, i);
		bool vr = is_valid(rhs, i);
		return (vl && vr) || (vl && lhs[i]) || (vr && rhs[i]);
	}
};

// Combine masks
template<Binop Op, Vec L, Vec R>
constexpr auto mask_join(const L lhs, const R rhs) noexcept
{
	if constexpr ( Op == And )
		return vec_valid_and<L,R>{lhs, rhs};
	else if constexpr ( Op == Or )
		return vec_valid_or<L,R>{lhs, rhs};
	else
	{
		if constexpr ( Masked<L> && Masked<R> )
		{
			auto ml = lhs.get_mask();
			auto mr = rhs.get_mask();
			return vec_valid<decltype(ml),decltype(mr)>
			{
				.ml = ml,
				.mr = mr,
			};
		}
		else if constexpr ( Masked<L> )
			return lhs.get_mask();
		else if constexpr ( Masked<R> )
			return rhs.get_mask();
		else
			static_assert(Op != Op, "at least one operand must be masked");
	}
}

// Mask a vector
template<Num T = double, Vec V, Vec Mask>
constexpr auto mask(
	const V data,
	const Mask mask) noexcept
{
	assert(data.ssize() == mask.ssize());
	if constexpr ( Masked<V> )
	{
		auto _data = data.get_data();
		auto _mask = data.get_mask();
		auto newmask = vec_valid<Mask,decltype(_mask)>{mask, _mask};
		return vec_masked<decltype(_data),decltype(newmask),T>
		{
			.data = _data,
			.mask = newmask,
		};
	}
	else
	{
		return vec_masked<V,Mask,T>
		{
			.data = data,
			.mask = mask,
		};
	}
}

// Mask a vector to exclude NAs
template<Num T = double, Vec V>
constexpr auto mask(const V data) noexcept
{
	return mask(data, ufunc<NotNA>(data));
}

// Vector subscripted at the given indices
template<Vec V, Vec Index, Num T = double>
struct vec_indexed
{
	V data;
	Index index;

	ptrdiff_t ssize() const noexcept
	{
		return index.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		if ( is_valid(index, i) )
			return coerce_cast<T>(data[index[i]]);
		else
			return na_value<T>();
	}
};

// Gather vector elements at given indices
template<Num T = double, Vec Index, Vec V>
constexpr auto gather(
	const Index index,
	const V data) noexcept
{
	if constexpr ( Masked<V> )
	{
		auto _data = data.get_data();
		auto _mask = data.get_mask();
		auto newdata = vec_indexed<decltype(_data),Index,T>
		{
			.data = _data,
			.index = index,
		};
		auto newmask = vec_indexed<decltype(_mask),Index,T>
		{
			.data = _mask,
			.index = index,
		};
		return mask<T>(newdata, newmask);
	}
	else
	{
		return vec_indexed<V,Index,T>
		{
			.data = data, 
			.index = index,
		};
	}
}

// Vector with elementwise unary transformation
template<Vec V, UnaryOp Tform, Num T = double>
struct vec_unop
{
	V data;
	Tform op;

	ptrdiff_t ssize() const noexcept
	{
		return data.ssize();
	}

	T operator[](ptrdiff_t i) const noexcept
	{
		return coerce_cast<T>(op(data[i]));
	}
};

// Transform with elementwise unary functor
template<Unop Op, Num T = double, Vec V, UnaryOp Tform = unop<Op,T>>
constexpr auto transform(
	const V data, 
	const Tform op = Tform{}) noexcept
{
	if constexpr ( Masked<V> )
	{
		auto _data = data.get_data();
		auto newdata = vec_unop<decltype(_data),Tform,T>
		{
			.data = _data,
			.op = op,
		};
		return mask<T>(newdata, data.get_mask());
	}
	else
	{
		return vec_unop<V,Tform,T>
		{
			.data = data, 
			.op = op,
		};
	}
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
	const Tform op = Tform{}) noexcept
{
	assert(lhs.ssize() == rhs.ssize());
	if constexpr ( Masked<L> && Masked<R> )
	{
		auto _lhs = lhs.get_data();
		auto _rhs = rhs.get_data();
		auto newdata = vec_binop<decltype(_lhs),decltype(_rhs),Tform,T>
		{
			.lhs = _lhs,
			.rhs = _rhs,
			.op = op,
		};
		auto newmask = mask_join<Op>(lhs, rhs);
		return mask<T>(newdata, newmask);
	}
	else if constexpr ( Masked<L> )
	{
		auto _lhs = lhs.get_data();
		auto newdata = vec_binop<decltype(_lhs),R,Tform,T>
		{
			.lhs = _lhs,
			.rhs = rhs,
			.op = op,
		};
		return mask<T>(newdata, lhs.get_mask());
	}
	else if constexpr ( Masked<R> )
	{
		auto _rhs = rhs.get_data();
		auto newdata = vec_binop<L,decltype(_rhs),Tform,T>
		{
			.lhs = lhs,
			.rhs = _rhs,
			.op = op,
		};
		return mask<T>(newdata, rhs.get_mask());
	}
	else
	{
		return vec_binop<L,R,Tform,T>
		{
			.lhs = lhs, 
			.rhs = rhs, 
			.op = op,
		};
	}
}

// Reduce vector elements with binary functor
template<Binop Op, Num T = double, Vec V, BinaryOp Reduce = binop<Op,T>>
T reduce(
	const V data,
	const Reduce op = Reduce{},
	const T init = Reduce::identity()) noexcept
{
	T accum = init;
	for ( ptrdiff_t i = 0; i < data.ssize(); ++i )
	{
		if ( data.is_valid(i) )
			accum = op(accum, coerce_cast<T>(data[i]));
	}
	return accum;
}

//// Vector operators
//---------------------
// Deferred math and arithmetic

// Universal unary functions for Vecs
template<Unop Op, Num T, Vec V>
constexpr Vec auto ufunc(V x) noexcept {
	return transform<Op,T>(x);
}

// Universal binary functions for Vecs
template<Binop Op, Num T, Vec L, Vec R>
constexpr Vec auto ufunc(L lhs, R rhs) noexcept {
	return transform<Op,T>(lhs, rhs);
}

// Vec + Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator+(L lhs, R rhs) noexcept {
	return ufunc<Add,T>(lhs, rhs);
}

// Vec - Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator-(L lhs, R rhs) noexcept {
	return ufunc<Sub,T>(lhs, rhs);
}

// Vec * Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator*(L lhs, R rhs) noexcept {
	return ufunc<Mul,T>(lhs, rhs);
}

// Vec / Vec
template<Num T = double, Vec L, Vec R>
constexpr Vec auto operator/(L lhs, R rhs) noexcept {
	return ufunc<Div,T>(lhs, rhs);
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
		assert(!is_null());
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const noexcept
	{
		assert(!is_null());
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	vec<T>& swap(const ptrdiff_t i, const ptrdiff_t j) noexcept
	{
		T xi = (*this)[i];
		(*this)[i] = (*this)[j];
		(*this)[j] = xi;
		return (*this);
	}

	vec<T>& fill(const T value) noexcept
	{
		return this->assign(rep<T>{value, len});
	}

	vec<T>& seqfill(const T start, const T step = 1) noexcept
	{
		return this->assign(seq<T>{start, len, step});
	}

	// Elementwise assignment
	template<Vec V>
	vec<T>& assign(const V src) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			if constexpr ( Masked<V> )
			{
				if ( !src.is_valid(i) )
					continue;
			}
			(*this)[i] = coerce_cast<T>(src[i]);
		}
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
		{
			if constexpr ( Masked<V> )
			{
				if ( !src.is_valid(i) )
					continue;
			}
			(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		}
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
			if constexpr ( Masked<V> )
			{
				if ( !src.is_valid(ii) )
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
			if constexpr ( Masked<V> )
			{
				if ( !src.is_valid(i) )
					continue;
			}
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
		return this->transform<Sub>(src);
	}

	// vec<T> *= Vec
	template<Vec V>
	vec<T>& operator*=(const V src) noexcept
	{
		return this->transform<Mul>(src);
	}

	// vec<T> /= Vec
	template<Vec V>
	vec<T>& operator/=(const V src) noexcept
	{
		return this->transform<Div>(src);
	}

	// Return a sliced view from b.start to b.stop
	vec<T> slice(const bounds b) const noexcept
	{
		assert(b.start <= b.stop);
		assert(0 <= b.start && b.start < len);
		assert(0 <= b.stop && b.stop <= len);
		return {
			.ptr = ptr + (stride * b.start), 
			.len = b.stop - b.start,
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

	ptrdiff_t ssize() const noexcept
	{
		return nrows * ncols;
	}

	T& operator[](const loc<2> i) noexcept
	{
		assert(0 <= i[0] && i[0] < nrows);
		assert(0 <= i[1] && i[1] < ncols);
		return ptr[row_stride * i[0] + col_stride * i[1]];
	}

	const T& operator[](const loc<2> i) const noexcept
	{
		assert(0 <= i[0] && i[0] < nrows);
		assert(0 <= i[1] && i[1] < ncols);
		return ptr[row_stride * i[0] + col_stride * i[1]];
	}

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
			.nrows = b.stop - b.start,
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
			.ncols = b.stop - b.start,
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

//// Ragged arrays
//-----------------
// Vectors of different lenghts

// A non-owning ragged array
template<Num T, Num Offset>
struct rag
{
	vec<T> x;
	vec<Offset> offset;

	ptrdiff_t ssize() const noexcept
	{
		return offset.len - 1;
	}

	vec<T> operator[](ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i + 1 < offset.len);
		return x.slice({offset[i], offset[i + 1]});
	}

	#ifdef USING_R
	static rag<T,Offset> from(SEXP x, SEXP offset) noexcept
	{
		return {
			.x = vec<T>::from(x),
			.offset = vec<Offset>::from(offset),
		};
	}
	#endif // USING_R

};

#endif // CARDINAL_CORE_CORE
