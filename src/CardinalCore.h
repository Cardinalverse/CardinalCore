#ifndef CARDINAL_CORE
#define CARDINAL_CORE

#include "R.h"
#include "Rinternals.h"
#include "search.h"

extern "C" {

// Search and selection
SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative_diff);

} // extern "C"

#endif // CARDINAL_CORE
