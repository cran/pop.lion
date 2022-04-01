#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP C_montecarlo(SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP,
                         SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_montecarlo", (DL_FUNC) &C_montecarlo, 14},
    {NULL, NULL, 0}
};

void R_init_pop_lion(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
