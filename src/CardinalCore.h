#ifndef CARDINAL_CORE
#define CARDINAL_CORE

#include "R.h"
#include "Rinternals.h"
#include "search.h"

extern "C" {

// Search and selection
//----------------------
SEXP do_sdiff(SEXP x, SEXP ref, SEXP relative);

} // extern "C"

#endif // CARDINAL_CORE
