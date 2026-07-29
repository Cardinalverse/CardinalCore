
#include "CardinalCore.h"

extern "C" {

//// Quicksort and Quickselect
//----------------------------

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP result = PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(x); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				REAL(result)[i] = diff(
					INTEGER_ELT(x, i),
					INTEGER_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			case REALSXP:
				REAL(result)[i] = diff(
					REAL_ELT(x, i),
					REAL_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			default:
				Rf_error("'x' and 'ref' must be integer or double");
		}
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qorder(SEXP x)
{
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			quick_order(
				as_vec<int>(result).fill(0, 1),
				as_vec<int>(x));
			break;
		case REALSXP:
			quick_order(
				as_vec<int>(result).fill(0, 1),
				as_vec<double>(x));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qselect(SEXP x, SEXP k)
{
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), LENGTH(k)));
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(k); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				INTEGER(result)[i] = quick_select(
					as_vec<int>(index).fill(0, 1),
					as_vec<int>(x),
					INTEGER_ELT(k, i));
				break;
			case REALSXP:
				REAL(result)[i] = quick_select(
					as_vec<int>(index).fill(0, 1),
					as_vec<double>(x),
					INTEGER_ELT(k, i));
				break;
			default:
				Rf_error("'x' must be integer or double");
		}
	}
	UNPROTECT(2);
	return result;
}

SEXP do_qmedian(SEXP x)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_median(as_vec<int>(x)));
			break;
		case REALSXP:
			return Rf_ScalarReal(quick_median(as_vec<double>(x)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_qmad(SEXP x, SEXP center, SEXP constant)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_mad(
				as_vec<int>(x),
				Rf_asReal(center),
				Rf_asReal(constant)));
		case REALSXP:
			return Rf_ScalarReal(quick_mad(
				as_vec<double>(x),
				Rf_asReal(center),
				Rf_asReal(constant)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_bsearch(
	SEXP query, 
	SEXP x, 
	SEXP tolerance, 
	SEXP relative, 
	SEXP nearest, 
	SEXP nomatch)
{
	if ( TYPEOF(query) != TYPEOF(x) )
		Rf_error("'query' and 'x' must have the same data type");
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(query)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			binary_search(
				INTEGER(result),
				as_vec<int>(query),
				as_vec<int>(x),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			binary_search(
				INTEGER(result),
				as_vec<double>(query),
				as_vec<double>(x),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

//// Matrix statistics
//---------------------

SEXP do_col_sums(SEXP x, SEXP num_threads)
{
	SEXP result = PROTECT(Rf_allocVector(REALSXP, Rf_ncols(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			chunk_apply(
				col_sums<int>{
					as_vec<double>(result).fill(),
					as_mat<int>(x)},
				Rf_asInteger(num_threads),
				XLENGTH(result));
			break;
		}
		case REALSXP:
		{
			chunk_apply(
				col_sums<double>{
					as_vec<double>(result).fill(),
					as_mat<double>(x)},
				Rf_asInteger(num_threads),
				XLENGTH(result));
			break;
		}
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
