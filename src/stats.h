#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include "core.h"
#include "kernels.h"

//// Stats
//---------
// Traits supporting streaming statistics

enum Summary {
	Sum,
	Prod,
	Maximum,
	Minimum,
	Mean,
	Var,
};

// Binary operation for simple reduction summaries
template<Summary S>
constexpr Binop summary_op() noexcept
{
	if constexpr ( S == Sum )
		return Add;
	else if constexpr ( S == Prod )
		return Mul;
	else if constexpr ( S == Maximum )
		return Max;
	else if constexpr ( S == Minimum )
		return Min;
	else
		static_assert(dependent_false_v<S>, "unsupported summary op");
}

// A Stat supports a streaming statistic with online updates
// - SHOULD also support s.update(Num x) -> Stat
// - SHOULD also support s.merge(Stat s) -> Stat
template<class S>
concept Stat = 
	requires (std::remove_cvref_t<S>& s)
	{
		{ s.get() } -> Num;
		{ s.nobs() } -> Num;
	};

// A Stats Vec supports vector access to streaming statistics
// - SHOULD also support s.get(ptrdiff_t i) -> Stat
// - SHOULD also support s.set(ptrdiff_t i, Stat s)
template<class S>
concept Stats = Vec<S> &&
	requires (std::remove_cvref_t<S>& s, ptrdiff_t i)
	{
		{ s.get_stats() } -> Vec;
		{ s.get_nobs() } -> Vec;
		{ s.nobs(i) } -> Num;
	};

//// Stats Vec operations
//-----------------------
// Merge and update streaming statistics

// Used to implement dst.merge(src)
template<Stats Dst, Stat Val>
Dst fill_stats(Dst dst, const Val value) noexcept
{
	for ( ptrdiff_t i = 0; i < dst.ssize(); ++i )
		dst.set(i, value);
	return dst;
}

// Used to implement dst.update(src)
template<Stats Dst, Vec Src>
Dst update_stats(Dst dst, const Src src) noexcept
{
	assert(dst.ssize() == src.ssize());
	for ( ptrdiff_t i = 0; i < dst.ssize(); ++i )
	{
		if ( is_valid(src, i) )
			dst.set(i, dst.get(i).update(src[i]));
	}
	return dst;
}

// Used to implement dst.merge(src)
template<Stats Dst, Stats Src>
Dst merge_stats(Dst dst, const Src src) noexcept
{
	assert(dst.ssize() == src.ssize());
	for ( ptrdiff_t i = 0; i < dst.ssize(); ++i )
		dst.set(i, dst.get(i).merge(src.get(i)));
	return dst;
}

// Used to implement dst.scatter_update(index, src)
template<Stats Dst, Vec Index, Vec Src>
Dst scatter_update_stats(Dst dst, const Index index, const Src src) noexcept
{
	assert(dst.ssize() == src.ssize());
	for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
	{
		if ( !is_valid(index, i) || !is_valid(src, i) )
			continue;
		auto ii = index[i];
		dst.set(ii, dst.get(ii).update(src[i]));
	}
	return dst;
}

// Used to implement dst.scatter_merge(index, src)
template<Stats Dst, Vec Index, Stats Src>
Dst scatter_merge_stats(Dst dst, const Index index, const Src src) noexcept
{
	assert(dst.ssize() == src.ssize());
	for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
	{
		if ( !is_valid(index, i) )
			continue;
		auto ii = index[i];
		dst.set(ii, dst.get(ii).merge(src.get(i)));
	}
	return dst;
}

//// Scalar statistics
//--------------------
// Streaming scalar statistics

// Simple reduction
template<Summary S, Num T, Num N>
struct stream_stat
{
	T stat = binop<summary_op<S>(),T>::identity();
	N n = 0;

	T get() const noexcept { return stat; }
	N nobs() const noexcept { return n; }

	stream_stat<S,T,N> update(const T x) const noexcept
	{
		return {
			.stat = ufunc<summary_op<S>(),T>(stat, x),
			.n = n + 1,
		};
	}

	stream_stat<S,T,N> merge(const stream_stat<S,T,N> s) const noexcept
	{
		return {
			.stat = ufunc<summary_op<S>(),T>(stat, s.stat),
			.n = n + s.n,
		};
	}

	template<Vec V>
	stream_stat<S,T,N> reduce(const V src) const noexcept
	{
		auto dst = *this;
		for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
		{
			if ( is_valid(src, i) )
				dst = dst.update(src[i]);
		}
		return dst;
	}
};

// Mean specialization
template<Num T, Num N>
struct stream_stat<Mean,T,N>
{
	T mean = na_value<T>();
	N n = 0;

	T get() const noexcept { return mean; }
	N nobs() const noexcept { return n; }

	stream_stat<Mean,T,N> update(const T x) const noexcept
	{
		if ( n > 0 ) {
			return {
				.mean = ((n * mean) + x) / (n + 1),
				.n = n + 1,
			};
		}
		else {
			return {
				.mean = x,
				.n = 1,
			};
		}
	}

	stream_stat<Mean,T,N> merge(const stream_stat<Mean,T,N> s) const noexcept
	{
		if ( n > 0 ) {
			if ( s.n > 0 )
			{
				return {
					.mean = ((n * mean) + (s.n * s.mean)) / (n + s.n),
					.n = n + s.n,
				};
			}
			else
				return (*this);
		}
		else {
			if ( s.n > 0 )
				return s;
			else
				return (*this);
		}
	}

	template<Vec V>
	stream_stat<Mean,T,N> reduce(const V src) const noexcept
	{
		T sout = n > 0 ? n * mean : 0.0;
		N nout = n;
		for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
		{
			if ( is_valid(src, i) ) {
				sout += src[i];
				++nout;
			}
		}
		return {
			.mean = sout / nout,
			.n = nout,
		};
	}
};

// Variance specialization
template<Num T, Num N>
struct stream_stat<Var,T,N>
{
	T var = na_value<T>();
	T mean = na_value<T>();
	N n = 0;

	T get() const noexcept { return var; }
	N nobs() const noexcept { return n; }

	// Online variance update from Welford (1962)
	stream_stat<Var,T,N> update(const T x) const noexcept
	{
		if ( n >= 2 ) {
			T m0 = mean;
			T m1 = ((n * mean) + x) / (n + 1);
			T ss0 = (n - 1) * var;
			T ss1 = ss0 + (x - m0) * (x - m1);
			return {
				.var = ss1 / n,
				.mean = m1,
				.n = n + 1,
			};
		}
		else if ( n >= 1 ) {
			return {
				.var = ((mean - x) * (mean - x)) / 2,
				.mean = ((n * mean) + x) / (n + 1),
				.n = n + 1,
			};
		}
		else {
			return {
				.var = var,
				.mean = x,
				.n = 1,
			};
		}
	}

	// Batch variance update from Chan, Golub, & LeVeque (1979)
	stream_stat<Var,T,N> merge(const stream_stat<Var,T,N> s) const noexcept
	{
		if ( n > 1 ) {
			if ( s.n > 1 )
			{
				T nn = static_cast<T>(n * s.n) / static_cast<T>(n + s.n);
				T mm = (mean - s.mean) * (mean - s.mean);
				T ssa = (n - 1) * var;
				T ssb = (s.n - 1) * s.var;
				T ssnew = ssa + ssb + (mm * nn);
				return {
					.var = ssnew / (n + s.n - 1),
					.mean = ((n * mean) + (s.n * s.mean)) / (n + s.n),
					.n = static_cast<N>(n + s.n),
				};
			}
			else if ( s.n == 1 )
				return this->update(s.mean);
			else
				return (*this);
		}
		else if ( n == 1 ) {
			if ( s.n > 1 )
				return s.update(mean);
			else if ( s.n == 1)
				return this->update(s.mean);
			else
				return (*this);
		}
		else {
			if ( s.n > 0 )
				return s;
			else
				return (*this);
		}
	}

	template<Vec V>
	stream_stat<Var,T,N> reduce(const V src) const noexcept
	{
		auto init = stream_stat<Mean,T,N>{mean, n};
		auto xbar = init.reduce(src);
		T ssout = 0;
		for ( ptrdiff_t i = 0; i < src.ssize(); ++i )
		{
			if ( is_valid(src, i) ) {
				T dx = src[i] - xbar.mean;
				ssout += (dx * dx);
			}
		}
		return {
			.var = ssout / (xbar.n - 1),
			.mean = xbar.mean,
			.n = xbar.n,
		};
	}
};

//// Vector statistics
//--------------------
// Streaming vector statistics

// Vector with streaming summary stats
// - Pairs reductions and numbers of observations
// - MUST have n.ssize() == stats.ssize()
template<Summary S, Num T, Num N>
struct stream_stats
{
	vec<T> stats{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 0 ? stats[i] : binop<summary_op<S>(),T>::identity();
	}

	vec<T> get_stats() const noexcept { return stats; }
	vec<N> get_nobs() const noexcept { return n; }
	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	stream_stat<S,T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.stat = stats[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, stream_stat<S,T,N> s) noexcept
	{
		stats[i] = s.stat;
		n[i] = s.n;
	}

	stream_stats<S,T,N> fill(stream_stat<S,T,N> value = {}) noexcept {
		return fill_stats(*this, value);
	}

	template<Vec V>
	stream_stats<S,T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_stats<S,T,N> merge(const stream_stats<S,T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index, Vec V>
	stream_stats<S,T,N> scatter_update(
		const Index index, 
		const V src) noexcept 
	{
		return scatter_update_stats(*this, index, src);
	}

	template<Vec Index>
	stream_stats<S,T,N> scatter_merge(
		const Index index, 
		const stream_stats<S,T,N> s) noexcept 
	{
		return scatter_merge_stats(*this, index, s);
	}

	template<Vec V>
	void reduce_into(const ptrdiff_t i, V x) noexcept {
		set(i, stream_stat<S,T,N>{}.reduce(x));
	}

	stream_stats<S,T,N> slice(bounds b) const noexcept
	{
		return {
			.stats = stats.slice(b),
			.n = n.slice(b),
		};
	}

	#ifdef USING_R
	static stream_stats<S,T,N> from(SEXP obj) noexcept
	{
		return {
			.stats = r_vec<T>(obj),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};


// Vector with streaming means
// - Pairs means and numbers of observations
// - MUST have n.ssize() == means.ssize()
template<Num T, Num N>
struct stream_stats<Mean,T,N>
{
	vec<T> means{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 0 ? means[i] : na_value<T>();
	}

	vec<T> get_stats() const noexcept { return means; }
	vec<N> get_nobs() const noexcept { return n; }
	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	stream_stat<Mean,T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.mean = means[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, stream_stat<Mean,T,N> s) noexcept
	{
		means[i] = s.mean;
		n[i] = s.n;
	}

	stream_stats<Mean,T,N> fill(stream_stat<Mean,T,N> value = {}) noexcept {
		return fill_stats(*this, value);
	}

	template<Vec V>
	stream_stats<Mean,T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_stats<Mean,T,N> merge(const stream_stats<Mean,T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index, Vec V>
	stream_stats<Mean,T,N> scatter_update(
		const Index index, 
		const V src) noexcept 
	{
		return scatter_update_stats(*this, index, src);
	}

	template<Vec Index>
	stream_stats<Mean,T,N> scatter_merge(
		const Index index, 
		const stream_stats<Mean,T,N> s) noexcept 
	{
		return scatter_merge_stats(*this, index, s);
	}

	template<Vec V>
	void reduce_into(const ptrdiff_t i, V x) noexcept {
		set(i, stream_stat<Mean,T,N>{}.reduce(x));
	}

	stream_stats<Mean,T,N> slice(bounds b) const noexcept
	{
		return {
			.means = means.slice(b),
			.n = n.slice(b),
		};
	}

	#ifdef USING_R
	static stream_stats<Mean,T,N> from(SEXP obj) noexcept
	{
		return {
			.means = r_vec<T>(obj),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

// Vector with streaming variances
// - Pairs variances and numbers of observations
// - MUST have n.ssize() == vars.ssize() == means.ssize()
template<Num T, Num N>
struct stream_stats<Var,T,N>
{
	vec<T> vars{};
	vec<T> means{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 1 ? vars[i] : na_value<T>();
	}

	vec<T> get_stats() const noexcept { return vars; }
	vec<N> get_nobs() const noexcept { return n; }
	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	stream_stat<Var,T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.var = vars[i],
			.mean = means[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, const stream_stat<Var,T,N> s) noexcept
	{
		vars[i] = s.var;
		means[i] = s.mean;
		n[i] = s.n;
	}

	stream_stats<Var,T,N> fill(stream_stat<Var,T,N> value = {}) noexcept {
		return fill_stats(*this, value);
	}

	template<Vec V>
	stream_stats<Var,T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_stats<Var,T,N> merge(const stream_stats<Var,T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index, Vec V>
	stream_stats<Var,T,N> scatter_update(
		const Index index, 
		const V src) noexcept 
	{
		return scatter_update_stats(*this, index, src);
	}

	template<Vec Index>
	stream_stats<Var,T,N> scatter_merge(
		const Index index, 
		const stream_stats<Var,T,N> s) noexcept 
	{
		return scatter_merge_stats(*this, index, s);
	}

	template<Vec V>
	void reduce_into(const ptrdiff_t i, V x) noexcept {
		set(i, stream_stat<Var,T,N>{}.reduce(x));
	}

	stream_stats<Var,T,N> slice(bounds b) const noexcept
	{
		return {
			.vars = vars.slice(b),
			.means = means.slice(b),
			.n = n.slice(b),
		};
	}

	#ifdef USING_R
	static stream_stats<Var,T,N> from(SEXP obj) noexcept
	{
		return {
			.vars = r_vec<T>(obj),
			.means = r_vec<T>(Rf_getAttrib(obj, Rf_install("mean"))),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

//// Multichannel statistics
//--------------------------
// Multiple streams of statistics

// Vectors with streaming summary stats
// - Consists of ssize() streams with nchannels() features
// - Where get(i) -> stream_stats with nchannels() elements
// - Where channel(j) -> stream_stats with ssize() elements
template<Summary S, Num T, Num N, Dim Channels = Rows>
struct nstream_stats
{
	mat<T> stats{};
	mat<N> n{};

	ptrdiff_t ssize() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.ncols();
		else
			return n.nrows();
	}

	ptrdiff_t nchannels() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.nrows();
		else
			return n.ncols();
	}

	stream_stats<S,T,N> channel(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {stats.row(i), n.row(i)};
		else
			return {stats.col(i), n.col(i)};
	}

	nstream_stats<S,T,N> fill(stream_stat<S,T,N> value = {}) const noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			get(i).fill(value);
		return *this;
	}

	stream_stats<S,T,N> get(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {stats.col(i), n.col(i)};
		else
			return {stats.row(i), n.row(i)};
	}

	nstream_stats<S,T,N,Channels> slice(bounds b) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {stats.slice_cols(b), n.slice_cols(b)};
		else
			return {stats.slice_rows(b), n.slice_rows(b)};
	}

	#ifdef USING_R
	static nstream_stats<S,T,N,Channels> from(SEXP obj) noexcept
	{
		return {
			.stats = r_mat<T>(obj),
			.n = r_mat<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

// Vectors with streaming means
// - Consists of ssize() streams with nchannels() features
// - Where get(i) -> stream_stats with nchannels() elements
// - Where channel(j) -> stream_stats with ssize() elements
template<Num T, Num N, Dim Channels>
struct nstream_stats<Mean,T,N,Channels>
{
	mat<T> means{};
	mat<N> n{};

	ptrdiff_t ssize() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.ncols();
		else
			return n.nrows();
	}

	ptrdiff_t nchannels() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.nrows();
		else
			return n.ncols();
	}

	stream_stats<Mean,T,N> channel(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {means.row(i), n.row(i)};
		else
			return {means.col(i), n.col(i)};
	}

	nstream_stats<Mean,T,N> fill(stream_stat<Mean,T,N> value = {}) const noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			get(i).fill(value);
		return *this;
	}

	stream_stats<Mean,T,N> get(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {means.col(i), n.col(i)};
		else
			return {means.row(i), n.row(i)};
	}

	nstream_stats<Mean,T,N,Channels> slice(bounds b) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {means.slice_cols(b), n.slice_cols(b)};
		else
			return {means.slice_rows(b), n.slice_rows(b)};
	}

	#ifdef USING_R
	static nstream_stats<Mean,T,N,Channels> from(SEXP obj) noexcept
	{
		return {
			.means = r_mat<T>(obj),
			.n = r_mat<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

// Vectors with streaming variances
// - Consists of ssize() streams with nchannels() features
// - Where get(i) -> stream_stats with nchannels() elements
// - Where channel(j) -> stream_stats with ssize() elements
template<Num T, Num N, Dim Channels>
struct nstream_stats<Var,T,N,Channels>
{
	mat<T> vars{};
	mat<T> means{};
	mat<N> n{};

	ptrdiff_t ssize() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.ncols();
		else
			return n.nrows();
	}

	ptrdiff_t nchannels() const noexcept
	{
		if constexpr ( Channels == Rows )
			return n.nrows();
		else
			return n.ncols();
	}

	stream_stats<Var,T,N> channel(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {vars.row(i), means.row(i), n.row(i)};
		else
			return {vars.col(i), means.col(i), n.col(i)};
	}

	nstream_stats<Var,T,N> fill(stream_stat<Var,T,N> value = {}) const noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			get(i).fill(value);
		return *this;
	}

	stream_stats<Var,T,N> get(ptrdiff_t i) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {vars.col(i), means.col(i), n.col(i)};
		else
			return {vars.row(i), means.row(i), n.row(i)};
	}

	nstream_stats<Var,T,N,Channels> slice(bounds b) const noexcept
	{
		if constexpr ( Channels == Rows )
			return {vars.slice_cols(b), means.slice_cols(b), n.slice_cols(b)};
		else
			return {vars.slice_rows(b), means.slice_rows(b), n.slice_rows(b)};
	}

	#ifdef USING_R
	static nstream_stats<Var,T,N,Channels> from(SEXP obj) noexcept
	{
		return {
			.vars = r_mat<T>(obj),
			.means = r_mat<T>(Rf_getAttrib(obj, Rf_install("mean"))),
			.n = r_mat<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

//// Aggregate
//------------
// Summarize a vector

// Aggregate a summary statistic at indices
template<Summary S, Num T, Num N, Vec V>
struct aggregate
{
	stream_stat<S,T,N> * stat; // out
	V x;                       // in

	void operator()(ptrdiff_t i) noexcept
	{
		if ( is_valid(x, i) )
			*stat = stat->update(coerce_cast<T>(x[i]));
	}

	void operator()(bounds b) noexcept
	{
		for ( ptrdiff_t i = b.start; i < b.stop; ++i )
			(*this)(i);
	}
};

//// Stats kernels
//----------------
// Summarize matrices

template<Summary S, Num T, Num N, Mat M>
struct col_stats
{
	stream_stats<S,T,N> stats; // out
	M x;                       // in

	ptrdiff_t ssize() const { return stats.ssize(); }

	void operator()()
	{
		if ( x.prefer_rows() )
			for ( ptrdiff_t i = 0; i < x.nrows(); ++i )
				stats.update(mask(x.row(i)));
		else
			for ( ptrdiff_t j = 0; j < x.ncols(); ++j )
				stats.reduce_into(j, mask(x.col(j)));
	}

	void operator()(bounds b, task ctx)
	{
		col_stats{stats.slice(b), x.slice_cols(b)}();
	}
};

template<Summary S, Num T, Num N, Mat M>
struct row_stats
{
	stream_stats<S,T,N> stats; // out
	M x;                       // in

	ptrdiff_t ssize() const { return stats.ssize(); }

	void operator()()
	{
		if ( x.prefer_rows() )
			for ( ptrdiff_t i = 0; i < x.nrows(); ++i )
				stats.reduce_into(i, mask(x.row(i)));
		else
			for ( ptrdiff_t j = 0; j < x.ncols(); ++j )
				stats.update(mask(x.col(j)));
	}

	void operator()(bounds b, task ctx)
	{
		row_stats{stats.slice(b), x.slice_rows(b)}();
	}
};

template<Summary S, Num T, Num N, Dim Channels, Mat M, Vec Index>
struct col_scatter_stats
{
	nstream_stats<S,T,N,Channels> stats; // out
	Index index;                         // in
	M x;                                 // in

	ptrdiff_t ssize() const { return stats.ssize(); }

	void operator()()
	{
		if ( x.prefer_rows() )
		{
			for ( ptrdiff_t i = 0; i < x.nrows(); ++i ) {
				if ( is_valid(index, i) )
					stats.channel(index[i]).update(mask(x.row(i)));
			}
		}
		else
			for ( ptrdiff_t j = 0; j < x.ncols(); ++j )
				stats.get(j).scatter_update(index, mask(x.col(j)));
	}

	void operator()(bounds b, task ctx)
	{
		col_scatter_stats{stats.slice(b), index, x.slice_cols(b)}();
	}
};

template<Summary S, Num T, Num N, Dim Channels, Mat M, Vec Index>
struct row_scatter_stats
{
	nstream_stats<S,T,N,Channels> stats; // out
	Index index;                         // in
	M x;                                 // in

	ptrdiff_t ssize() const { return stats.ssize(); }

	void operator()()
	{
		if ( x.prefer_rows() )
			for ( ptrdiff_t i = 0; i < x.nrows(); ++i )
				stats.get(i).scatter_update(index, mask(x.row(i)));
		else
		{
			for ( ptrdiff_t j = 0; j < x.ncols(); ++j ) {
				if ( is_valid(index, j) )
					stats.channel(index[j]).update(mask(x.col(j)));
			}
		}
	}

	void operator()(bounds b, task ctx)
	{
		row_scatter_stats{stats.slice(b), index, x.slice_rows(b)}();
	}
};


#endif // CARDINAL_CORE_STATS
