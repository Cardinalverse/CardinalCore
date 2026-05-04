#ifndef CARDINAL_CORE_PRELUDE
#define CARDINAL_CORE_PRELUDE

//// Structs
//-----------

template<typename T>
struct matrix 
{
	T * data;
	size_t nrows;
	size_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;

	inline T at(size_t i, size_t j)
	{
		return data[row_stride * i + col_stride * j];
	}

	inline const T * row(size_t i)
	{
		return data + (row_stride * i);
	}

	inline const T * col(size_t j)
	{
		return data + (col_stride * j);
	}
};

template<typename T>
struct scaling
{
	T * center;
	T * scale;
	int * group;
	size_t ngroups;
	ptrdiff_t group_stride;

	inline T center_at(size_t i, size_t grp)
	{
		if ( ngroups <= 1 )
			return center[i];
		else
			return center[i + (group_stride * group[grp])];
	}

	inline T scale_at(size_t i, size_t grp)
	{
		if ( ngroups <= 1 )
			return scale[i];
		else
			return scale[i + (group_stride * group[grp])];
	}
};

//// Comparison
//--------------

inline bool isIncomparable(int x)
{
	return x == NA_INTEGER;
}

inline bool isIncomparable(double x)
{
	return ISNAN(x);
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
template<typename T> inline
double diff(T x, T ref, bool relative = false)
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

#define LESSER(x, y) (diff((x), (y)) < 0)
#define GREATER(x, y) (diff((x), (y)) > 0)
#define LESSER_EQUAL(x, y) (diff((x), (y)) <= 0)
#define GREATER_EQUAL(x, y) (diff((x), (y)) >= 0)
#define EQUAL(x, y) (diff((x), (y)) == 0)
#define NOT_EQUAL(x, y) (diff((x), (y)) != 0)

//// Utility
//-----------

template<typename T>
void fill_buffer(
	T * buffer,      // buffer to fill
	size_t size,     // size of buffer
	T init = 0,      // value to use to initialize
	T increment = 0, // increment init for each item
	int stride = 1)
{
	for ( size_t i = 0; i < size; i++ )
	{
		buffer[stride * i] = init;
		init += increment;
	}
}

#endif // CARDINAL_CORE_PRELUDE
