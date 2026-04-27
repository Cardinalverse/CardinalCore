
#include "CardinalCore.h"
#include "search.h"

extern "C" {

SEXP do_sdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
	{
		if ( LENGTH(x) > 1 && LENGTH(y) > 1 )
			Rf_error("'x' and 'ref' must be equal length or length 1")
	}
	double result = sdiff(
		Rf_asReal(x),
		Rf_asReal(ref),
		Rf_asLogical(relative));
	return Rf_ScalarReal(result);
}

} // extern "C"
