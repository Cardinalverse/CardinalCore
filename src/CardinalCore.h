#ifndef CARDINAL_CORE
#define CARDINAL_CORE

// #ifdef NDEBUG
// #undef NDEBUG
// #endif

#define R_NO_REMAP
#include "R.h"
#include "Rinternals.h"
#include "order.h"
#include "search.h"
#include "kernels.h"

extern "C" {

//// Sort and order
//------------------
SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative);
SEXP do_qorder(SEXP x);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(SEXP x, SEXP center, SEXP constant);

//// Search and nearest neighbors
//--------------------------------
SEXP do_bsearch(SEXP x, SEXP data, SEXP tolerance, 
	SEXP relative, SEXP nearest, SEXP nomatch);

//// Matrix statistics
//---------------------
SEXP do_col_sums(SEXP x, SEXP num_threads);

//// Test expressions
//--------------------
SEXP do_test_expression(SEXP x, SEXP index);

} // extern "C"

#endif // CARDINAL_CORE
