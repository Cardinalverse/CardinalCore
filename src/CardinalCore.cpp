
#include "CardinalCore.h"
#include "search.h"

extern "C" {

SEXP do_qdiff(
	SEXP x,
	SEXP ref,
	SEXP relative_diff)
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
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
