
#include "CardinalCore.h"
#include "search.h"

extern "C" {

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP result = PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			do_qdiff_impl<int>(
				INTEGER_RO(x),
				INTEGER_RO(ref),
				REAL(result),
				LENGTH(x),
				Rf_asLogical(relative));
			break;
		case REALSXP:
			do_qdiff_impl<double>(
				REAL_RO(x),
				REAL_RO(ref),
				REAL(result),
				LENGTH(x),
				Rf_asLogical(relative));
			break;
		default:
			Rf_error("'x' and 'ref' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qselect(SEXP x, SEXP k)
{
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), LENGTH(k)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			do_qselect_impl<int,int,int>(
				INTEGER_RO(x),
				LENGTH(x),
				INTEGER_RO(k),
				LENGTH(k),
				INTEGER(result));
			break;
		case REALSXP:
			do_qselect_impl<double,int,int>(
				REAL_RO(x),
				LENGTH(x),
				INTEGER_RO(k),
				LENGTH(k),
				REAL(result));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qsort(SEXP x)
{
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			quick_sort<int,int>(
				INTEGER_RO(x),
				0,
				LENGTH(x),
				INTEGER(result),
				true);
			break;
		case REALSXP:
			quick_sort<double,int>(
				REAL_RO(x),
				0,
				LENGTH(x),
				INTEGER(result),
				true);
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qmedian(SEXP x)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_median<int,int>(
				INTEGER_RO(x),
				LENGTH(x)));
			break;
		case REALSXP:
			return Rf_ScalarReal(quick_median<double,int>(
				REAL_RO(x),
				LENGTH(x)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_qmad(SEXP x, SEXP center, SEXP scale)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_mad<int,int>(
				INTEGER_RO(x),
				LENGTH(x),
				Rf_asReal(center),
				Rf_asReal(scale)));
		case REALSXP:
			return Rf_ScalarReal(quick_mad<double,int>(
				REAL_RO(x),
				LENGTH(x),
				Rf_asReal(center),
				Rf_asReal(scale)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_bsearch(SEXP x, SEXP data, SEXP tolerance, 
	SEXP relative, SEXP nearest, SEXP nomatch)
{
	if ( TYPEOF(x) != TYPEOF(data) )
		Rf_error("'x' and 'data' must have the same data type");
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			do_bsearch_impl<int,int>(
				INTEGER_RO(x),
				LENGTH(x),
				INTEGER_RO(data),
				LENGTH(data),
				INTEGER(result),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			do_bsearch_impl<double,int>(
				REAL_RO(x),
				LENGTH(x),
				REAL_RO(data),
				LENGTH(data),
				INTEGER(result),
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

} // extern "C"
