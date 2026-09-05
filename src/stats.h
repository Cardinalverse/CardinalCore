#ifndef CARDINAL_CORE_STATS
#define CARDINAL_CORE_STATS

#include "core.h"

//// Stats
//---------
// Traits supporting streaming statistics

// A Stats Vec supports vector access to streaming statistics
template<class S>
concept Stats = Vec<S> &&
	requires (std::remove_cvref_t<S>& s, ptrdiff_t i)
	{
		{ s.nobs(i) } -> Num;
	};

//// Stats operations
//-------------------
// Merge and update streaming statistics

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

// Used to implement dst.scatter(index, src)
template<Stats Dst, Vec Index, Stats Src>
Dst scatter_stats(Dst dst, const Index index, const Src src) noexcept
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

// Summary (online updates)
template<Binop Op, Num T, Num N>
struct stream_stat
{
	T stat = binop<Op,T>::identity();
	N n = 0;

	stream_stat<Op,T,N> update(const T x) const noexcept
	{
		return {
			.stat = ufunc<Op,T>(stat, x),
			.n = n + 1,
		};
	}

	stream_stat<Op,T,N> merge(const stream_stat<Op,T,N> s) const noexcept
	{
		return {
			.stat = ufunc<Op,T>(stat, s.stat),
			.n = n + s.n,
		};
	}
};

// Mean (online updates)
template<Num T, Num N>
struct stream_mean
{
	T mean = na_value<T>();
	N n = 0;

	stream_mean<T,N> update(const T x) const noexcept
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

	stream_mean<T,N> merge(const stream_mean<T,N> s) const noexcept
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
};

// Variance (online updates)
template<Num T, Num N>
struct stream_var
{
	T var = na_value<T>();
	T mean = na_value<T>();
	N n = 0;

	// Online variance update from Welford (1962)
	stream_var<T,N> update(const T x) const noexcept
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
	stream_var<T,N> merge(const stream_var<T,N> s) const noexcept
	{
		if ( n > 1 ) {
			if ( s.n > 1 )
			{
				N nn = (n * s.n) / (n + s.n);
				T mm = (mean - s.mean) * (mean - s.mean);
				T ssa = (n - 1) * var;
				T ssb = (s.n - 1) * s.var;
				T ssnew = ssa + ssb + (mm * nn);
				return {
					.var = ssnew / (n + s.n - 1),
					.mean = ((n * mean) + (s.n * s.mean)) / (n + s.n),
					.n = n + s.n,
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
};

//// Vector statistics
//--------------------
// Streaming vector statistics

// A vector with streaming summary stats
// - Pairs reductions and numbers of observations
// - MUST have n.ssize() == stats.ssize()
template<Binop Op, Num T, Num N>
struct stream_stats
{
	vec<T> stats{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 0 ? stats[i] : binop<Op,T>::identity();
	}

	stream_stat<Op,T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.stat = stats[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, stream_stat<Op,T,N> s) noexcept
	{
		stats[i] = s.stat;
		n[i] = s.n;
	}

	stream_stats<Op,T,N>& fill(stream_stat<Op,T,N> value = {}) noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			set(i, value);
		return (*this);
	}

	template<Vec V>
	stream_stats<Op,T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_stats<Op,T,N> merge(const stream_stats<Op,T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index>
	stream_stats<Op,T,N> scatter(
		const Index index, 
		const stream_stats<Op,T,N> s) noexcept 
	{
		return scatter_stats(*this, index, s);
	}

	#ifdef USING_R
	static stream_stats<Op,T,N> from(SEXP obj) noexcept
	{
		return {
			.stats = r_vec<T>(obj),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};


// A vector with streaming means
// - Pairs means and numbers of observations
// - MUST have n.ssize() == means.ssize()
template<Num T, Num N>
struct stream_means
{
	vec<T> means{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 0 ? means[i] : na_value<T>();
	}

	stream_mean<T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.mean = means[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, stream_mean<T,N> s) noexcept
	{
		means[i] = s.mean;
		n[i] = s.n;
	}

	stream_means<T,N>& fill(stream_mean<T,N> value = {}) noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			set(i, value);
		return (*this);
	}

	template<Vec V>
	stream_means<T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_means<T,N> merge(const stream_means<T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index>
	stream_means<T,N> scatter(
		const Index index, 
		const stream_means<T,N> s) noexcept 
	{
		return scatter_stats(*this, index, s);
	}

	#ifdef USING_R
	static stream_means<T,N> from(SEXP obj) noexcept
	{
		return {
			.means = r_vec<T>(obj),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

// A vector with streaming variance
// - Pairs variances and numbers of observations
// - MUST have n.ssize() == vars.ssize() == means.ssize()
template<Num T, Num N>
struct stream_vars
{
	vec<T> vars{};
	vec<T> means{};
	vec<N> n{};

	ptrdiff_t ssize() const noexcept { return n.ssize(); }

	N nobs(ptrdiff_t i) const noexcept { return n[i]; }

	T operator[](ptrdiff_t i) const noexcept
	{
		return n[i] > 1 ? vars[i] : na_value<T>();
	}

	stream_var<T,N> get(const ptrdiff_t i) const noexcept
	{
		return {
			.var = vars[i],
			.mean = means[i],
			.n = n[i],
		};
	}

	void set(const ptrdiff_t i, const stream_var<T,N> s) noexcept
	{
		vars[i] = s.var;
		means[i] = s.mean;
		n[i] = s.n;
	}

	stream_vars<T,N>& fill(stream_var<T,N> value = {}) noexcept
	{
		for ( ptrdiff_t i = 0; i < ssize(); ++i )
			set(i, value);
		return (*this);
	}

	template<Vec V>
	stream_vars<T,N> update(const V x) noexcept {
		return update_stats(*this, x);
	}

	stream_vars<T,N> merge(const stream_vars<T,N> s) noexcept {
		return merge_stats(*this, s);
	}

	template<Vec Index>
	stream_vars<T,N> scatter(
		const Index index, 
		const stream_vars<T,N> s) noexcept 
	{
		return scatter_stats(*this, index, s);
	}

	#ifdef USING_R
	static stream_vars<T,N> from(SEXP obj) noexcept
	{
		return {
			.vars = r_vec<T>(obj),
			.means = r_vec<T>(Rf_getAttrib(obj, Rf_install("mean"))),
			.n = r_vec<N>(Rf_getAttrib(obj, Rf_install("nobs"))),
		};
	}
	#endif // USING_R
};

#endif // CARDINAL_CORE_STATS
