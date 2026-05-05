#ifndef CARDINAL_CORE_PRELUDE
#define CARDINAL_CORE_PRELUDE

//// Data Ptr (mutable)
//---------------------

template<typename T, typename Object>
T * data_ptr(Object x);

template<> inline
int * data_ptr<int,SEXP>(SEXP x)
{
	return INTEGER(x);
}

template<> inline
double * data_ptr<double,SEXP>(SEXP x)
{
	return REAL(x);
}

//// Data Ptr (immutable)
//-----------------------

template<typename T, typename Object>
const T * data_ptr_const(Object x);

template<> inline
const int * data_ptr_const<int,SEXP>(SEXP x)
{
	return INTEGER_RO(x);
}

template<> inline
const double * data_ptr_const<double,SEXP>(SEXP x)
{
	return REAL_RO(x);
}

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

	matrix(SEXP x)
	{
		if ( x != R_NilValue )
		{
			data = data_ptr<T,SEXP>(x);
			nrows = Rf_nrows(x);
			ncols = Rf_ncols(x);
			row_stride = 1;
			col_stride = Rf_nrows(x);
		}
		else
		{
			data = nullptr;
			nrows = 0;
			ncols = 0;
			row_stride = 0;
			col_stride = 0;
		}
	}

	inline bool valid() const
	{
		return data != nullptr;
	}

	inline T at(size_t i, size_t j) const
	{
		return data[row_stride * i + col_stride * j];
	}

	inline const T * row(size_t i) const
	{
		return data + (row_stride * i);
	}

	inline const T * col(size_t j) const
	{
		return data + (col_stride * j);
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
