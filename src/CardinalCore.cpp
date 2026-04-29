
#include "CardinalCore.h"
#include "search.h"

extern "C" {

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative_diff)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( XLENGTH(x) != XLENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP result = PROTECT(Rf_allocVector(REALSXP, XLENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			do_qdiff_impl<int>(
				INTEGER_RO(x),
				INTEGER_RO(ref),
				REAL(result),
				static_cast<size_t>(XLENGTH(x)),
				static_cast<bool>(Rf_asLogical(relative_diff)));
			break;
		case REALSXP:
			do_qdiff_impl<double>(
				REAL_RO(x),
				REAL_RO(ref),
				REAL(result),
				static_cast<size_t>(XLENGTH(x)),
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
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(k)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			do_qselect_impl<int,int,int>(
				INTEGER_RO(x),
				0,
				XLENGTH(x),
				INTEGER_RO(k),
				XLENGTH(k),
				INTEGER(result));
			break;
		case REALSXP:
			do_qselect_impl<double,int,int>(
				REAL_RO(x),
				0,
				XLENGTH(x),
				INTEGER_RO(k),
				XLENGTH(k),
				REAL(result));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
