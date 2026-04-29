
#include <R_ext/Rdynload.h>

#include "CardinalCore.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef callMethods[] = {
	// sort and search
	CALLDEF(do_qdiff, 3),
	CALLDEF(do_qselect, 2),
	CALLDEF(do_qsort, 1),
	CALLDEF(do_qmedian, 1),
	CALLDEF(do_qmad, 3),
	{NULL, NULL, 0}
};

void R_init_CardinalCore(DllInfo * info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

} // extern "C"
