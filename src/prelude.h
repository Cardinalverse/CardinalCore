#ifndef CARDINAL_CORE_PRELUDE
#define CARDINAL_CORE_PRELUDE

//// Data Ptr (mutable)
//---------------------
// Get a mutable data pointer from an object

template<typename T>
T * data_ptr(SEXP x);

template<> inline
int * data_ptr<int>(SEXP x)
{
	return INTEGER(x);
}

template<> inline
double * data_ptr<double>(SEXP x)
{
	return REAL(x);
}

//// Data Ptr (immutable)
//-----------------------
// Get an immutable data pointer from an object

template<typename T>
const T * data_ptr_const(SEXP x);

template<> inline
const int * data_ptr_const<int>(SEXP x)
{
	return INTEGER_RO(x);
}

template<> inline
const double * data_ptr_const<double>(SEXP x)
{
	return REAL_RO(x);
}

//// Structs
//-----------
// Containers for common object types

struct slice 
{
	ptrdiff_t begin;
	ptrdiff_t end;

	inline size_t size() const
	{
		if ( end > begin )
			return static_cast<size_t>(end - begin);
		else
			return 0;
	}
};

template<typename T>
struct vctr 
{
	T * data;
	size_t len;
	ptrdiff_t stride;

	inline T at(const ptrdiff_t i) const
	{
		return data[stride * i];
	}

	inline slice all_elements() const
	{
		return {0, static_cast<size_t>(len)};
	}

};

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
			data = data_ptr<T>(x);
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

	inline T at(
		const ptrdiff_t i, 
		const ptrdiff_t j) const
	{
		return data[(row_stride * i) + (col_stride * j)];
	}

	inline vctr<T> row_vctr(const ptrdiff_t i) const
	{
		return {data + (row_stride * i), ncols, col_stride};
	}

	inline vctr<T> col_vctr(const ptrdiff_t i) const
	{
		return {data + (col_stride * i), nrows, row_stride};
	}

	inline slice all_rows() const
	{
		return {0, static_cast<ptrdiff_t>(nrows)};
	}

	inline slice all_cols() const
	{
		return {0, static_cast<ptrdiff_t>(ncols)};
	}
};

template<typename T>
matrix<T> as_matrix(SEXP x)
{
	if ( x != R_NilValue )
	{
		return {
			.data = data_ptr<T>(x),
			.nrows = Rf_nrows(x),
			.ncols = Rf_ncols(x),
			.row_stride = 1,
			.col_stride = Rf_nrows(x)
		};
	}
	else
	{
		return {
			.data = nullptr,
			.nrows = 0,
			.ncols = 0,
			.row_stride = 0,
			.col_stride = 0
		};
	}
}

//// Incomparables
//-----------------
// Handle incomparable values (NAs and NaNs)

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

//// Comparison
//--------------
// compute signed absolute or relative difference
// * safe to use with incomparables (NAs and NaNs)
// * incomparables sort last/highest (NA >> Inf)
// * incomparables sort equal to each other (NA == NA)
// returns: the difference
template<typename T> inline
double diff(
	const T x, 
	const T ref, 
	const bool relative = false)
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
	const size_t size,     // size of buffer
	T init = 0,      // value to use to initialize
	const T increment = 0, // increment init for each item
	const int stride = 1)
{
	for ( size_t i = 0; i < size; i++ )
	{
		buffer[stride * i] = init;
		init += increment;
	}
}

#endif // CARDINAL_CORE_PRELUDE
