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

//// Min/Max
//------------

#define MIN2(x, y) ((x) < (y) ? (x) : (y))
#define MIN3(x, y, z) (MIN2(MIN2((x), (y)), (z)))
#define MAX2(x, y) ((x) > (y) ? (x) : (y))
#define MAX3(x, y, z) (MAX2(MAX2((x), (y)), (z)))

//// Structs
//-----------
// Containers for common object types

enum Axis {
	Rows,
	Columns
};

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

struct groups
{
	int * index;
	int ngroups;
};

template<typename T>
struct vctr 
{
	T * ptr;
	size_t len;
	ptrdiff_t stride;

	inline T at(const ptrdiff_t i) const
	{
		return ptr[stride * i];
	}

	inline slice all_elements() const
	{
		return {0, static_cast<ptrdiff_t>(len)};
	}

};

template<typename T>
struct matrix 
{
	T * ptr;
	size_t nrows;
	size_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;

	inline T at(
		const ptrdiff_t i, 
		const ptrdiff_t j) const
	{
		return ptr[(row_stride * i) + (col_stride * j)];
	}

	inline size_t dim(Axis axis) const
	{
		switch(axis) {
			case Rows:
				return nrows;
			case Columns:
				return ncols;
			default:
				return 0;
		}
	}

	inline vctr<T> row_vctr(const ptrdiff_t i) const
	{
		return {ptr + (row_stride * i), ncols, col_stride};
	}

	inline vctr<T> col_vctr(const ptrdiff_t i) const
	{
		return {ptr + (col_stride * i), nrows, row_stride};
	}

	inline slice all_rows() const
	{
		return {0, static_cast<ptrdiff_t>(nrows)};
	}

	inline slice all_cols() const
	{
		return {0, static_cast<ptrdiff_t>(ncols)};
	}

	inline slice all_along(Axis axis) const
	{
		return {0, static_cast<ptrdiff_t>(dim(axis))};
	}

};

struct chunks
{
	size_t i;
	size_t len;
	int nchunks;

	inline slice next()
	{
		// compute remaining items
		size_t n = len - i;
		// exit early if no remaining chunks or items
		if ( nchunks == 0 || n == 0 )
		{
			return {
				.begin = static_cast<ptrdiff_t>(len),
				.end = static_cast<ptrdiff_t>(len)
			};
		}
		// estimate size of next chunk
		size_t chunksize = n / nchunks;
		// adjust chunksize for remaining items
		if ( chunksize < n * nchunks )
			++chunksize;
		else if ( chunksize > n )
			chunksize = n;
		// create slice
		slice chunk = {
			.begin = static_cast<ptrdiff_t>(i),
			.end = static_cast<ptrdiff_t>(i + chunksize)
		};
		// update members
		i += chunksize;
		--nchunks;
		// return slice
		return chunk;
	}
};

//// R objects
//--------------
// Initialize struct from R object

template<typename T>
vctr<T> as_vctr(SEXP x)
{
	if ( x != R_NilValue )
	{
		return {
			.ptr = data_ptr<T>(x),
			.len = static_cast<size_t>(XLENGTH(x)),
			.stride = 1
		};
	}
	else
	{
		return {
			.ptr = nullptr,
			.len = 0,
			.stride = 0
		};
	}
}

template<typename T>
matrix<T> as_matrix(SEXP x)
{
	if ( x != R_NilValue )
	{
		return {
			.ptr = data_ptr<T>(x),
			.nrows = static_cast<size_t>(Rf_nrows(x)),
			.ncols = static_cast<size_t>(Rf_ncols(x)),
			.row_stride = 1,
			.col_stride = Rf_nrows(x)
		};
	}
	else
	{
		return {
			.ptr = nullptr,
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

//// Infinities
//-----------------
// Define infinite values
#define POS_INF R_PosInf
#define NEG_INF R_NegInf 

//// Comparison
//--------------
// Comparisons handling incomparables (NAs and NaNs)

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
		return POS_INF; // NAs sort last so (x - ref) => +Inf
	else if ( isIncomparable(ref) )
		return NEG_INF; // NAs sort last so (x - ref) => -Inf
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
	T * buffer,            // buffer to fill
	const size_t size,     // size of buffer
	T init = 0,            // value to use to initialize
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
