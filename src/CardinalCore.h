#ifndef CARDINAL_CORE
#define CARDINAL_CORE

#include "R.h"
#include "Rinternals.h"
#include "search.h"

extern "C" {

// Search and selection
SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qsort(SEXP x);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(SEXP x, SEXP center, SEXP constant);
SEXP do_bsearch(SEXP x, SEXP data, SEXP tolerance, 
	SEXP relative, SEXP nearest, SEXP nomatch);

} // extern "C"

#endif // CARDINAL_CORE
