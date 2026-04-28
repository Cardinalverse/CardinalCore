#ifndef CARDINAL_CORE_SEARCH
#define CARDINAL_CORE_SEARCH

//// Comparison
//--------------

inline bool isIncomparable(int x)
{
	return x == NA_INTEGER;
}

inline bool isIncomparable(double x)
{
	return ISNA(x) || ISNAN(x);
}

// signed absolute or relative difference
template<typename T>
double qdiff(T x, T ref, bool relative_diff = false)
{
	if ( isIncomparable(x) && isIncomparable(ref) )
		return 0.0;
	else if ( isIncomparable(x) )
		return R_PosInf; // NAs sort last so (x - ref) => +Inf
	else if ( isIncomparable(ref) )
		return R_NegInf; // NAs sort last so (x - ref) => -Inf
	else
	{
		if ( relative_diff )
			return static_cast<double>(x - ref) / ref;
		else
			return static_cast<double>(x - ref);
	}
}

template<typename T>
void do_qdiff_impl(
	const T * x,
	const T * ref,
	double * out_result,
	size_t n,
	bool relative_diff)
{
	for ( size_t i = 0; i < n; ++i )
	{
		out_result[i] = qdiff<T>(x[i], ref[i], relative_diff);
	}
}

#endif // CARDINAL_CORE_SEARCH
