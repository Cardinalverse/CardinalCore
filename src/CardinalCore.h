#ifndef CARDINAL_CORE
#define CARDINAL_CORE

#define R_NO_REMAP
#include "R.h"
#include "Rinternals.h"
#include "search.h"
#include "stats.h"

extern "C" {

// Search and selection
SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qorder(SEXP x);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(SEXP x, SEXP center, SEXP constant);
SEXP do_bsearch(SEXP x, SEXP data, SEXP tolerance, 
	SEXP relative, SEXP nearest, SEXP nomatch);
SEXP do_col_sums(SEXP x, SEXP num_threads);
SEXP do_col_scatter_sums(SEXP x, SEXP group, SEXP ngroups, SEXP num_threads);

} // extern "C"

#endif // CARDINAL_CORE
