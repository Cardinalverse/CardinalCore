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

// Get type of a Vec's elements
template<class V>
using typeof_vec = std::remove_cvref_t<decltype(
	std::declval<const std::remove_cvref_t<V>&>()[std::declval<ptrdiff_t>()])>;

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
constexpr bool is_valid(V& v, ptrdiff_t i)
{
	if constexpr ( Masked<V> )
		return v.is_valid(i);
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

// Static assert false workaround
template<class T>
constexpr bool dependent_false = false;

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
		count += !is_valid(x, i) || is_na(x[i]);
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
	Identity,
	// Logic
	Not, IsNA, NotNA,
	// Math
	Abs, Sign, Log, Log2, Log1p, Exp, Exp2, Expm1,
};

template<Unop Op, Num T = double>
constexpr T ufunc(T x) noexcept
{
	// Identity
	if constexpr ( Op == Identity )
		return x;
	// Logic
	else if constexpr ( Op == Not )
		return !x;
	else if constexpr ( Op == IsNA )
		return is_na(x);
	else if constexpr ( Op == NotNA )
		return !is_na(x);
	// Math
	else if constexpr ( Op == Abs )
		return std::abs(x);
	else if constexpr ( Op == Sign )
		return (x > 0) - (x < 0);
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
		static_assert(dependent_false<T>, "unsupported unary op");
}

template<Unop Op, Num T = double>
struct unop
{
	T operator()(T x) const noexcept
	{
		return ufunc<Op,T>(x);
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
	Add, Sub, Mul, Div, Pow, Max, Min
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
		return lhs && rhs;
	else if constexpr ( Op == Or )
		return lhs || rhs;
	// Compare
	else if constexpr ( Op == Eq )
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
	else if constexpr ( Op == Pow )
		return std::pow(lhs, rhs);
	else if constexpr ( Op == Max )
		return lhs > rhs ? lhs : rhs;
	else if constexpr ( Op == Min )
		return lhs < rhs ? lhs : rhs;
	// Not implemented
	else
		static_assert(dependent_false<T>, "unsupported binary op");
}

template<Binop Op, Num T = double>
struct binop
{
	static T identity() noexcept
	{
		// Logic
		if constexpr ( Op == And )
			return true;
		else if constexpr ( Op == Or )
			return false;
		// Arithmetic
		else if constexpr ( Op == Add )
			return 0;
		else if constexpr ( Op == Mul )
			return 1;
		else if constexpr ( Op == Max )
		{
			if constexpr ( std::numeric_limits<T>::has_infinity )
				return -std::numeric_limits<T>::infinity();
			else
				return std::numeric_limits<T>::lowest();
		}
		else if constexpr ( Op == Min )
		{
			if constexpr ( std::numeric_limits<T>::has_infinity )
				return +std::numeric_limits<T>::infinity();
			else
				return std::numeric_limits<T>::max();
		}
		// Not implemented
		else
			static_assert(dependent_false<T>, "unsupported reduction");
	}

	T operator()(T lhs, T rhs) const noexcept
	{
		return ufunc<Op,T>(lhs, rhs);
	}
};

//// Vector expressions
//---------------------
// Lazy expressions on vectors

// Vector with elementwise unary transformation
template<Vec V, UnaryOp Tform, Num T = double>
struct vec_unop
{
	V x;
	Tform op;

	ptrdiff_t ssize() const noexcept { return x.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		if ( is_valid(x, i) )
			return coerce_cast<T>(op(x[i]));
		else
			return na_value<T>();
	}
};

// Vector with elementwise binary transformation
template<Vec L, Vec R, BinaryOp Tform, Num T = double>
struct vec_binop
{
	L lhs;
	R rhs;
	Tform op;

	ptrdiff_t ssize() const noexcept { return lhs.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		if ( is_valid(lhs, i) && is_valid(rhs, i) )
			return coerce_cast<T>(op(lhs[i], rhs[i]));
		else
			return na_value<T>();
	}
};

// Vector subscripted at the given indices
template<Vec V, Vec Index, Num T = double>
struct vec_indexed
{
	V data;
	Index index;

	ptrdiff_t ssize() const noexcept { return index.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		if ( is_valid(index, i) )
			return coerce_cast<T>(data[index[i]]);
		else
			return na_value<T>();
	}
};

// Vector with validity mask
template<Vec V, Vec Mask, Num T = double>
struct vec_masked
{
	V data;
	Mask mask;

	ptrdiff_t ssize() const noexcept { return data.ssize(); }

	constexpr V get_data() const noexcept { return data; }

	constexpr Mask get_mask() const noexcept { return mask; }

	T operator[](ptrdiff_t i) const noexcept
	{
		return is_valid(i) ? coerce_cast<T>(data[i]) : na_value<T>();
	}

	bool is_valid(ptrdiff_t i) const noexcept
	{
		return coerce_cast<bool>(mask[i]);
	}
};

// Mask with ternary AND logic
template<Vec L, Vec R, Num T = bool>
struct kleene_and
{
	L lhs;
	R rhs;

	ptrdiff_t ssize() const noexcept { return lhs.ssize(); }

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
struct kleene_or
{
	L lhs;
	R rhs;

	ptrdiff_t ssize() const noexcept { return lhs.ssize(); }

	// Valid if both valid or if either is valid-TRUE
	T operator[](ptrdiff_t i) const noexcept
	{
		bool vl = is_valid(lhs, i);
		bool vr = is_valid(rhs, i);
		return (vl && vr) || (vl && lhs[i]) || (vr && rhs[i]);
	}
};

// Combine masks with boolean logic
template<Vec L, Vec R>
constexpr Vec auto join_masks(const L ml, const R mr) noexcept
{
	auto op = binop<And,bool>{};
	return vec_binop<L,R,decltype(op),bool>
	{
		.lhs = ml,
		.rhs = mr,
		.op = op,
	};
}

// Combine masks from Masked Vecs with boolean logic
template<Vec L, Vec R>
constexpr Vec auto join_masks_of(const L lhs, const R rhs) noexcept
{
	if constexpr ( Masked<L> && Masked<R> )
		return join_masks(lhs.get_mask(), rhs.get_mask());
	else if constexpr ( Masked<L> )
		return lhs.get_mask();
	else if constexpr ( Masked<R> )
		return rhs.get_mask();
	else
		static_assert(dependent_false<L>, "neither operand is masked");
}

// Combine masks from Masked Vecs with ternary logic
template<Binop Op, Vec L, Vec R>
constexpr Vec auto join_masks_of(const L lhs, const R rhs) noexcept
{
	if constexpr ( Masked<L> || Masked<R> )
	{
		if constexpr ( Op == And )
			return kleene_and<L,R>{lhs, rhs};
		else if constexpr ( Op == Or )
			return kleene_or<L,R>{lhs, rhs};
		else
			return join_masks_of(lhs, rhs);
	}
	else
		static_assert(dependent_false<L>, "neither operand is masked");
}

// Mask a Vec
template<Vec V, Vec Mask>
constexpr Vec auto mask(const V data, const Mask mask) noexcept
{
	assert(data.ssize() == mask.ssize());
	if constexpr ( Masked<V> )
	{
		auto _data = data.get_data();
		auto _mask = data.get_mask();
		auto newmask = join_masks(mask, _mask);
		return vec_masked<decltype(_data),decltype(newmask),typeof_vec<V>>
		{
			.data = _data,
			.mask = newmask,
		};
	}
	else
	{
		return vec_masked<V,Mask,typeof_vec<V>>
		{
			.data = data,
			.mask = mask,
		};
	}
}

// Mask a vector to exclude NAs (and NaNs)
template<Vec V>
constexpr Vec auto mask(const V data) noexcept
{
	auto op = unop<NotNA,typeof_vec<V>>{};
	auto _mask = vec_unop<V,decltype(op),bool>{data, op};
	return mask(data, _mask);
}

// Unmask a Vec
template<Vec V>
constexpr Vec auto unmask(const V data) noexcept
{
	if constexpr ( Masked<V> )
		return data.get_data();
	else
		return data;
}

//// Vectors generators
//---------------------
// Lazy repetitions and sequences

// Generator vector repeating a constant
template<Num T = double>
struct rep
{
	T value;
	ptrdiff_t len;

	constexpr ptrdiff_t ssize() const noexcept { return len; }

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

	constexpr ptrdiff_t ssize() const noexcept { return len; }

	constexpr T operator[](const ptrdiff_t i) const noexcept
	{
		assert(0 <= i && i < len);
		return start + (i * step);
	}
};

//// Vector transformations
//-------------------------
// Gather, transform, and reduce

// Gather vector elements at given indices
template<Num T = double, Vec Index, Vec V>
constexpr Vec auto gather(const Index index, const V data) noexcept
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
		return mask(newdata, newmask);
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

// Transform with elementwise unary functor
template<Num T = double, Vec V, UnaryOp Tform>
constexpr Vec auto transform(const V x, const Tform op) noexcept
{
	if constexpr ( Masked<V> )
	{
		auto _x = x.get_data();
		auto out = vec_unop<decltype(_x),Tform,T>
		{
			.x = _x,
			.op = op,
		};
		return mask(out, x.get_mask());
	}
	else
	{
		return vec_unop<V,Tform,T>
		{
			.x = x, 
			.op = op,
		};
	}
}

// Transform with elementwise unop
template<Unop Op, Num T = double, Vec V>
constexpr Vec auto transform(const V x) noexcept {
	return transform<T>(x, unop<Op,T>{});
}

// Transform with elementwise binary functor
template<Num T = double, Vec L, Vec R, BinaryOp Tform>
constexpr Vec auto transform(const L lhs, const R rhs, const Tform op) noexcept
{
	assert(lhs.ssize() == rhs.ssize());
	if constexpr ( Masked<L> || Masked<R> )
	{
		auto _lhs = unmask(lhs);
		auto _rhs = unmask(rhs);
		auto out = vec_binop<decltype(_lhs),decltype(_rhs),Tform,T>
		{
			.lhs = _lhs,
			.rhs = _rhs,
			.op = op,
		};
		return mask(out, join_masks_of(lhs, rhs));
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

// Transform with elementwise binop
template<Binop Op, Num T = double, Vec L, Vec R>
constexpr Vec auto transform(const L lhs, const R rhs) noexcept
{
	if constexpr ( Masked<L> || Masked<R> )
		return mask(
			transform(unmask(lhs), unmask(rhs), binop<Op,T>{}),
			join_masks_of<Op>(lhs, rhs));
	else
		return transform(lhs, rhs, binop<Op,T>{});
}

// Transform with elementwise binop and scalar RHS
template<Binop Op, Num T = double, Vec L, Num R>
constexpr Vec auto transform(const L lhs, const R rhs) noexcept
{
	auto _rhs = rep<T>{coerce_cast<T>(rhs), lhs.ssize()};
	return transform<Op,T>(lhs, _rhs);
}

// Transform with elementwise binop and scalar LHS
template<Binop Op, Num T = double, Num L, Vec R>
constexpr Vec auto transform(const L lhs, const R rhs) noexcept
{
	auto _lhs = rep<T>{coerce_cast<T>(lhs), rhs.ssize()};
	return transform<Op,T>(_lhs, rhs);
}

// Reduce vector to a scalar with a binary functor
template<Num T = double, Vec V, BinaryOp Reduce>
T reduce(const V x, const Reduce op, const T init) noexcept
{
	T accum = init;
	for ( ptrdiff_t i = 0; i < x.ssize(); ++i )
	{
		if ( is_valid(x, i) )
			accum = op(accum, coerce_cast<T>(x[i]));
	}
	return accum;
}

// Reduce vector to a scalar with a binop
template<Binop Op, Num T = double, Vec V>
T reduce(const V x) noexcept {
	return reduce<T>(x, binop<Op,T>{}, binop<Op,T>::identity());
}

//// Vector unary ops
//--------------------
// Universal functions with one vector

// Math
template<Vec V>
constexpr Vec auto abs(V x) noexcept {
	return transform<Abs,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto sign(V x) noexcept {
	return transform<Sign,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto log(V x) noexcept {
	return transform<Log,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto log2(V x) noexcept {
	return transform<Log2,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto log1p(V x) noexcept {
	return transform<Log1p,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto exp(V x) noexcept {
	return transform<Exp,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto exp2(V x) noexcept {
	return transform<Exp2,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Vec auto expm1(V x) noexcept {
	return transform<Expm1,typeof_vec<V>>(x);
}

// Summary
template<Vec V>
constexpr Num auto sum(V x) noexcept {
	return reduce<Add,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Num auto prod(V x) noexcept {
	return reduce<Mul,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Num auto all(V x) noexcept {
	return reduce<And,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Num auto any(V x) noexcept {
	return reduce<Or,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Num auto min(V x) noexcept {
	return reduce<Min,typeof_vec<V>>(x);
}
template<Vec V>
constexpr Num auto max(V x) noexcept {
	return reduce<Max,typeof_vec<V>>(x);
}

//// Vector binary ops
//---------------------
// Universal functions with two vectors

// Powers
template<Vec L, Vec R>
constexpr Vec auto pow(L lhs, R rhs) noexcept {
	return transform<Pow,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto pow(L lhs, R rhs) noexcept {
	return transform<Pow,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto pow(L lhs, R rhs) noexcept {
	return transform<Pow,typeof_vec<R>>(lhs, rhs);
}

// Minima
template<Vec L, Vec R>
constexpr Vec auto pmin(L lhs, R rhs) noexcept {
	return transform<Min,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto pmin(L lhs, R rhs) noexcept {
	return transform<Min,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto pmin(L lhs, R rhs) noexcept {
	return transform<Min,typeof_vec<R>>(lhs, rhs);
}

// Maxima
template<Vec L, Vec R>
constexpr Vec auto pmax(L lhs, R rhs) noexcept {
	return transform<Max,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto pmax(L lhs, R rhs) noexcept {
	return transform<Max,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto pmax(L lhs, R rhs) noexcept {
	return transform<Max,typeof_vec<R>>(lhs, rhs);
}

// Operator+
template<Vec L, Vec R>
constexpr Vec auto operator+(L lhs, R rhs) noexcept {
	return transform<Add,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator+(L lhs, R rhs) noexcept {
	return transform<Add,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator+(L lhs, R rhs) noexcept {
	return transform<Add,typeof_vec<L>>(lhs, rhs);
}

// Operator-
template<Vec L, Vec R>
constexpr Vec auto operator-(L lhs, R rhs) noexcept {
	return transform<Sub,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator-(L lhs, R rhs) noexcept {
	return transform<Sub,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator-(L lhs, R rhs) noexcept {
	return transform<Sub,typeof_vec<R>>(lhs, rhs);
}

// Operator*
template<Vec L, Vec R>
constexpr Vec auto operator*(L lhs, R rhs) noexcept {
	return transform<Mul,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator*(L lhs, R rhs) noexcept {
	return transform<Mul,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator*(L lhs, R rhs) noexcept {
	return transform<Mul,typeof_vec<R>>(lhs, rhs);
}

// Operator/
template<Vec L, Vec R>
constexpr Vec auto operator/(L lhs, R rhs) noexcept {
	return transform<Div,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator/(L lhs, R rhs) noexcept {
	return transform<Div,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator/(L lhs, R rhs) noexcept {
	return transform<Div,typeof_vec<R>>(lhs, rhs);
}

// Operator&
template<Vec L, Vec R>
constexpr Vec auto operator&(L lhs, R rhs) noexcept {
	return transform<And,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator&(L lhs, R rhs) noexcept {
	return transform<And,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator&(L lhs, R rhs) noexcept {
	return transform<And,typeof_vec<R>>(lhs, rhs);
}

// Operator|
template<Vec L, Vec R>
constexpr Vec auto operator|(L lhs, R rhs) noexcept {
	return transform<Or,typeof_vec<L>>(lhs, rhs);
}
template<Vec L, Num R>
constexpr Vec auto operator|(L lhs, R rhs) noexcept {
	return transform<Or,typeof_vec<L>>(lhs, rhs);
}
template<Num L, Vec R>
constexpr Vec auto operator|(L lhs, R rhs) noexcept {
	return transform<Or,typeof_vec<R>>(lhs, rhs);
}

//// Vectors
//-----------
// 1D array operations

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner is responsible for data validity
template<Num T>
struct vec
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	ptrdiff_t ssize() const noexcept { return len; }

	T& operator[](const ptrdiff_t i) noexcept
	{
		assert(ptr != nullptr);
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const noexcept
	{
		assert(ptr != nullptr);
		assert(0 <= i && i < len);
		return ptr[stride * i];
	}

	// Compare elements at i and j
	T compare(const ptrdiff_t i, const ptrdiff_t j) const noexcept
	{
		T lhs = (*this)[i];
		T rhs = (*this)[j];
		if ( !is_na(lhs) && !is_na(rhs) )
			return lhs - rhs;
		else
			return is_na(lhs) - is_na(rhs);
	}

	// Swap items at i and j
	void swap(const ptrdiff_t i, const ptrdiff_t j) noexcept
	{
		T xi = (*this)[i];
		(*this)[i] = (*this)[j];
		(*this)[j] = xi;
	}

	// Elementwise assignment
	template<Vec V>
	vec<T>& assign(const V src) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			if ( is_valid(src, i) )
				(*this)[i] = coerce_cast<T>(src[i]);
		}
		return (*this);
	}

	// Fill with constant value
	vec<T>& fill(const T value) noexcept {
		return this->assign(rep<T>{value, len});
	}

	// Fill with NAs
	vec<T>& fill_na() noexcept {
		return this->assign(rep<T>{na_value<T>(), len});
	}

	// Fill with sequential values
	vec<T>& fill_seq() noexcept {
		return this->assign(seq<T>{0, len, 1});
	}

	// Elementwise in-place unary transformations
	template<UnaryOp Tform>
	vec<T>& transform(const Tform op) noexcept
	{
		for ( ptrdiff_t i = 0; i < len; ++i )
			(*this)[i] = op((*this)[i]);
		return (*this);
	}

	// Elementwise in-place unop
	template<Unop Op>
	vec<T>& transform() noexcept {
		return transform(unop<Op,T>{});
	}

	// Elementwise in-place binary transformations
	template<Vec V, BinaryOp Tform>
	vec<T>& transform(const V src, const Tform op) noexcept
	{
		assert(src.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			if ( is_valid(src, i) )
				(*this)[i] = op((*this)[i], coerce_cast<T>(src[i]));
		}
		return (*this);
	}

	// Elementwise in-place binop
	template<Binop Op, Vec V>
	vec<T>& transform(const V src) noexcept {
		return transform(src, binop<Op,T>{});
	}

	// Elementwise in-place binop with a scalar
	template<Binop Op, Num N>
	vec<T>& transform(const N src) noexcept {
		auto _src = rep<T>{coerce_cast<T>(src), len};
		return transform(_src, binop<Op,T>{});
	}

	// Assign (*this)[i] = src[index[i]] for i in index
	template<Binop Op = Rhs, Vec Index, Vec V>
	vec<T>& gather(const Index index, const V src) noexcept
	{
		assert(index.ssize() == len);
		for ( ptrdiff_t i = 0; i < len; ++i )
		{
			if ( !is_valid(index, i) )
				continue;
			auto ii = index[i];
			if ( is_valid(src, ii) )
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
			if ( !is_valid(index, i) )
				continue;
			auto ii = index[i];
			if ( is_valid(src, ii) )
				(*this)[ii] = ufunc<Op,T>((*this)[ii], coerce_cast<T>(src[i]));
		}
		return (*this);
	}

	// Operator+=
	template<class S>
	vec<T>& operator+=(const S src) noexcept {
		return this->transform<Add>(src);
	}
	
	// Operator-=
	template<class S>
	vec<T>& operator-=(const S src) noexcept {
		return this->transform<Sub>(src);
	}

	// Operator*=
	template<class S>
	vec<T>& operator*=(const S src) noexcept {
		return this->transform<Mul>(src);
	}

	// Operator/=
	template<class S>
	vec<T>& operator/=(const S src) noexcept {
		return this->transform<Div>(src);
	}

	// Operator&=
	template<class S>
	vec<T>& operator&=(const S src) noexcept {
		return this->transform<And>(src);
	}

	// Operator|=
	template<class S>
	vec<T>& operator|=(const S src) noexcept {
		return this->transform<Or>(src);
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

	ptrdiff_t ssize() const noexcept { return nrows * ncols; }

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
