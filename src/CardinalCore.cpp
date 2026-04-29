
#include "CardinalCore.h"
#include "search.h"

extern "C" {

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative_diff)
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
				static_cast<size_t>(LENGTH(x)),
				static_cast<bool>(Rf_asLogical(relative_diff)));
			break;
		case REALSXP:
			do_qdiff_impl<double>(
				REAL_RO(x),
				REAL_RO(ref),
				REAL(result),
				static_cast<size_t>(LENGTH(x)),
				static_cast<bool>(Rf_asLogical(relative_diff)));
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
				0,
				LENGTH(x),
				INTEGER_RO(k),
				LENGTH(k),
				INTEGER(result));
			break;
		case REALSXP:
			do_qselect_impl<double,int,int>(
				REAL_RO(x),
				0,
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

} // extern "C"
