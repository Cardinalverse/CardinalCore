#ifndef CARDINAL_CORE_COMPARE
#define CARDINAL_CORE_COMPARE

//// Comparison
//--------------

#define LESSER(x, y) (qdiff((x), (y)) < 0)
#define GREATER(x, y) (qdiff((x), (y)) > 0)
#define LESSER_EQUAL(x, y) (qdiff((x), (y)) <= 0)
#define GREATER_EQUAL(x, y) (qdiff((x), (y)) >= 0)
#define EQUAL(x, y) (qdiff((x), (y)) == 0)
#define NOT_EQUAL(x, y) (qdiff((x), (y)) != 0)

inline bool isIncomparable(int x)
{
	return x == NA_INTEGER;
}

inline bool isIncomparable(double x)
{
	return ISNA(x) || ISNAN(x);
}

template<typename T>
T mkIncomparable();

template<> inline
int mkIncomparable<int>()
{
	return NA_INTEGER;
}

template<> inline
double mkIncomparable<double>()
{
	return NA_REAL;
}

// compute signed absolute or relative difference
// * safe to use with incomparables (NAs and NaNs)
// * incomparables sort last/highest (NA >> Inf)
// * incomparables sort equal to each other (NA == NA)
// returns: the difference
template<typename T>
double qdiff(T x, T ref, bool relative = false)
{
	if ( isIncomparable(x) && isIncomparable(ref) )
		return 0.0;
	else if ( isIncomparable(x) )
		return R_PosInf; // NAs sort last so (x - ref) => +Inf
	else if ( isIncomparable(ref) )
		return R_NegInf; // NAs sort last so (x - ref) => -Inf
	else
	{
		if ( relative )
			return static_cast<double>(x - ref) / ref;
		else
			return static_cast<double>(x - ref);
	}
}

#endif // CARDINAL_CORE_COMPARE
