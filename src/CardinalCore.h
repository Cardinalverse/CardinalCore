#ifndef CARDINAL_CORE
#define CARDINAL_CORE

#include "R.h"
#include "Rinternals.h"
#include "search.h"

extern "C" {

// Search and selection
SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative_diff);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qsort(SEXP x);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(SEXP x, SEXP center, SEXP constant);

} // extern "C"

#endif // CARDINAL_CORE
