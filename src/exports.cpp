
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "order.h"
#include "search.h"
#include "kernels.h"

extern "C" {

//// Quicksort and Quickselect
//----------------------------

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP dx = PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(x); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				REAL(dx)[i] = diff(
					INTEGER_ELT(x, i),
					INTEGER_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			case REALSXP:
				REAL(dx)[i] = diff(
					REAL_ELT(x, i),
					REAL_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			default:
				Rf_error("'x' and 'ref' must be integer or double");
		}
	}
	UNPROTECT(1);
	return dx;
}

SEXP do_qorder(SEXP x)
{
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			qsort_index(
				vec<int>::from(index).seqfill(),
				vec<int>::from(x));
			break;
		case REALSXP:
			qsort_index(
				vec<int>::from(index).seqfill(),
				vec<double>::from(x));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return index;
}

SEXP do_qselect(SEXP x, SEXP k)
{
	SEXP order = PROTECT(Rf_allocVector(TYPEOF(x), LENGTH(k)));
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(k); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				INTEGER(order)[i] = qselect_index(
					vec<int>::from(index).seqfill(),
					vec<int>::from(x),
					INTEGER_ELT(k, i));
				break;
			case REALSXP:
				REAL(order)[i] = qselect_index(
					vec<int>::from(index).seqfill(),
					vec<double>::from(x),
					INTEGER_ELT(k, i));
				break;
			default:
				Rf_error("'x' must be integer or double");
		}
	}
	UNPROTECT(2);
	return order;
}

SEXP do_qmedian(SEXP x)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(qmedian(vec<int>::from(x)));
			break;
		case REALSXP:
			return Rf_ScalarReal(qmedian(vec<double>::from(x)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_qmad(SEXP x, SEXP center, SEXP constant)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(qmad(
				vec<int>::from(x),
				Rf_asReal(center),
				Rf_asReal(constant)));
		case REALSXP:
			return Rf_ScalarReal(qmad(
				vec<double>::from(x),
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
				vec<int>::from(query),
				vec<int>::from(x),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			binary_search(
				INTEGER(result),
				vec<double>::from(query),
				vec<double>::from(x),
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
			compute(
				col_sums<int>{
					vec<double>::from(result).fill(),
					mat<int>::from(x)},
				Rf_asInteger(num_threads));
			break;
		}
		case REALSXP:
		{
			compute(
				col_sums<double>{
					vec<double>::from(result).fill(),
					mat<double>::from(x)},
				Rf_asInteger(num_threads));
			break;
		}
	}
	UNPROTECT(1);
	return result;
}

//// Test expressions
//--------------------

SEXP do_test_expression(SEXP x, SEXP index)
{
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(index)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			test_expression(
				vec<int>::from(result),
				vec<int>::from(x),
				vec<int>::from(index));
			break;
		}
		case REALSXP:
		{
			test_expression(
				vec<double>::from(result),
				vec<double>::from(x),
				vec<int>::from(index));
			break;
		}
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
