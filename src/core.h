#ifndef CARDINAL_CORE_CORE
#define CARDINAL_CORE_CORE

#include <cassert>
#include <cstddef>
#include <cmath>

//// Memory and errors
//---------------------
// Allocation and errors handled by managing runtime

#define SAFE_ALLOC R_Calloc
#define SAFE_FREE R_Free
#define SAFE_ERROR Rf_error

//// Data Pointer
//----------------
// Get a mutable data pointer from a managed object

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

//// Structs
//-----------
// Containers for common object types

// Index bounds
// - The interval is half-open: [start, stop)
struct bounds 
{
	ptrdiff_t start;
	ptrdiff_t stop;

	inline ptrdiff_t len() const
	{
		return stop - start;
	}
};

// A non-owning strided vector
// - Owner is responsible for managing memory
// - Owner MUST guarantee len >= 0
// - Callers MAY choose to handle stride < 0
template<typename T>
struct vec 
{
	T * ptr;
	ptrdiff_t len;
	ptrdiff_t stride;

	T& operator[](const ptrdiff_t i)
	{
		return ptr[stride * i];
	}

	const T& operator[](const ptrdiff_t i) const
	{
		return ptr[stride * i];
	}

	bounds all_elements() const
	{
		return {0, len};
	}

	vec<T> subset(bounds b) const
	{
		return {
			.ptr = ptr + (stride * b.start), 
			.len = b.len(),
			.stride = stride,
		};
	}
};

// A non-owning strided matrix
// - Owner is responsible for managing memory
// - Owner MUST guarantee nrows >= 0 and ncols >= 0
// - Callers MAY choose to handle row_stride < 0 and col_stride < 0
template<typename T>
struct mat 
{
	T * ptr;
	ptrdiff_t nrows;
	ptrdiff_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;

	vec<T> row(const ptrdiff_t i) const
	{
		return {
			.ptr = ptr + (row_stride * i), 
			.len = ncols, 
			.stride = col_stride
		};
	}

	vec<T> col(const ptrdiff_t i) const
	{
		return {
			.ptr = ptr + (col_stride * i), 
			.len = nrows, 
			.stride = row_stride
		};
	}

	bounds all_rows() const
	{
		return {0, nrows};
	}

	bounds all_cols() const
	{
		return {0, ncols};
	}

	mat<T> subset_rows(const bounds b) const
	{
		return {
			.ptr = ptr + (row_stride * b.start),
			.nrows = b.len(),
			.ncols = ncols,
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}

	mat<T> subset_cols(const bounds b) const
	{
		return {
			.ptr = ptr + (col_stride * b.start),
			.nrows = nrows,
			.ncols = b.len(),
			.row_stride = row_stride,
			.col_stride = col_stride,
		};
	}
};

//// R objects
//--------------
// Initialize struct from R object

template<typename T>
vec<T> as_vec(SEXP x)
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

template<typename T>
mat<T> as_mat(SEXP x)
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

//// Incomparables
//-----------------
// Handle incomparable values (NAs and NaNs)

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

inline bool isIncomparable(int x)
{
	return x == NA_INTEGER;
}

inline bool isIncomparable(double x)
{
	return ISNAN(x);
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
// - safe to use with incomparables (NAs and NaNs)
// - incomparables sort last/highest (NA >> Inf)
// - incomparables sort equal to each other (NA == NA)
// returns: the difference
template<typename T>
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
// Common idioms

#define MIN2(x, y) ((x) < (y) ? (x) : (y))
#define MIN3(x, y, z) (MIN2(MIN2((x), (y)), (z)))
#define MAX2(x, y) ((x) > (y) ? (x) : (y))
#define MAX3(x, y, z) (MAX2(MAX2((x), (y)), (z)))

template<typename T>
void fill_buffer(
	T * buffer,            // buffer to fill
	const ptrdiff_t size,     // size of buffer
	const T init = 0,            // value to use to initialize
	const T increment = 0,
	const ptrdiff_t stride = 1) // increment init for each item
{
	for ( ptrdiff_t i = 0; i < size; i += stride )
	{
		buffer[i] = init + (i * increment);
	}
}

#endif // CARDINAL_CORE_CORE
