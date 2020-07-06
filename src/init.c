#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "mcCNV.h"

static const R_CallMethodDef callMethods[] = {
  {"calcMLCN", (DL_FUNC) &calcMLCN, 5},
  {NULL, NULL, 0}
};

void R_init_mcCNV(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}