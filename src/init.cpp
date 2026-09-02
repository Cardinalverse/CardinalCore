
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {

//// Sort and order
//------------------
SEXP do_qorder(SEXP x);
SEXP do_qselect(SEXP x, SEXP k);
SEXP do_qmedian(SEXP x);
SEXP do_qmad(SEXP x, SEXP center, SEXP constant);

//// Search and nearest neighbors
//--------------------------------
SEXP do_bsearch(
	SEXP query,
	SEXP table,
	SEXP tolerance,
	SEXP relative,
	SEXP referent,
	SEXP nomatch);
SEXP do_rsearch(
	SEXP query,
	SEXP table,
	SEXP tolerance,
	SEXP relative,
	SEXP referent,
	SEXP nomatch);
SEXP do_kdtree_build(SEXP table);
SEXP do_kdtree_range_search(
	SEXP query,
	SEXP tree,
	SEXP tolerance,
	SEXP relative,
	SEXP referent,
	SEXP num_threads);
SEXP do_kdtree_knn_search(
	SEXP query,
	SEXP tree,
	SEXP k,
	SEXP p,
	SEXP num_threads);

//// Signal processing
//---------------------
SEXP do_filt1_mean(SEXP y, SEXP k);
SEXP do_filt1_conv(SEXP y, SEXP w);

//// Peak processing
//------------------
SEXP do_peaks_find(SEXP y, SEXP k);
SEXP do_peaks_snrs(SEXP y, SEXP k, SEXP method, SEXP wlen);
SEXP do_peaks_prominences(SEXP y, SEXP k, SEXP wlen);
SEXP do_peaks_widths(SEXP y, SEXP x, SEXP k, SEXP fmax);
SEXP do_peaks_areas(SEXP y, SEXP x, SEXP k);
SEXP do_peaks_summary(
	SEXP y, 
	SEXP x, 
	SEXP k, 
	SEXP method, 
	SEXP wlen, 
	SEXP fmax);

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
	CALLDEF(do_qselect, 2),
	CALLDEF(do_qorder, 1),
	CALLDEF(do_qmedian, 1),
	CALLDEF(do_qmad, 3),
	// Search and nearest neighbors
	CALLDEF(do_bsearch, 6),
	CALLDEF(do_rsearch, 6),
	CALLDEF(do_kdtree_build, 1),
	CALLDEF(do_kdtree_range_search, 6),
	CALLDEF(do_kdtree_knn_search, 5),
	// Signal processing
	CALLDEF(do_filt1_mean, 2),
	CALLDEF(do_filt1_conv, 2),
	// Peak processing
	CALLDEF(do_peaks_find, 2),
	CALLDEF(do_peaks_snrs, 4),
	CALLDEF(do_peaks_prominences, 3),
	CALLDEF(do_peaks_widths, 4),
	CALLDEF(do_peaks_areas, 3),
	CALLDEF(do_peaks_summary, 6),
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
