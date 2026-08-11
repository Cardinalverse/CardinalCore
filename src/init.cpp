
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {

//// Sort and order
//------------------
SEXP do_qdiff(
	SEXP x,
	SEXP ref,
	SEXP relative);
SEXP do_qorder(SEXP x);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(
	SEXP x,
	SEXP center,
	SEXP constant);

//// Search and nearest neighbors
//--------------------------------
SEXP do_bsearch(
	SEXP x,
	SEXP data,
	SEXP tolerance,
	SEXP relative,
	SEXP nearest,
	SEXP nomatch);
SEXP do_kdtree_build(SEXP data);
SEXP do_kdtree_range_search(
	SEXP queries,
	SEXP tree,
	SEXP tolerance,
	SEXP relative,
	SEXP num_threads);

//// Matrix statistics
//---------------------
SEXP do_col_sums(SEXP x, SEXP num_threads);

//// Test expressions
//--------------------
SEXP do_test_expression(SEXP x, SEXP index);

//// Register with R
//--------------------
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] = {
	// Sort and order
	CALLDEF(do_qdiff, 3),
	CALLDEF(do_qselect, 2),
	CALLDEF(do_qorder, 1),
	CALLDEF(do_qmedian, 1),
	CALLDEF(do_qmad, 3),
	// Search and nearest neighbors
	CALLDEF(do_bsearch, 6),
	CALLDEF(do_kdtree_build, 1),
	CALLDEF(do_kdtree_range_search, 5),
	// Matrix statistics
	CALLDEF(do_col_sums, 2),
	// Test expressions
	CALLDEF(do_test_expression, 2),
	{NULL, NULL, 0}
};

void R_init_CardinalCore(DllInfo * info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

} // extern "C"
