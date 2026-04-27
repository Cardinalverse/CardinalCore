
#include <R_ext/Rdynload.h>

#include "CardinalCore.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

static const R_CallMethodDef callMethods[] = {
	// sort and search
	CALLDEF(do_sdiff, 3),
	{NULL, NULL, 0}
};

void R_init_CardinalCore(DllInfo * info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

} // extern "C"
