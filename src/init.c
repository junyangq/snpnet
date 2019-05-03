#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP BEDMatrixPlus__extract_matrix(SEXP, SEXP, SEXP);
extern SEXP BEDMatrixPlus__extract_columns(SEXP, SEXP);
extern SEXP BEDMatrixPlus__new(SEXP, SEXP, SEXP);
extern SEXP BEDMatrixPlus__multiply_residuals(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callEntries[] = {
    {"BEDMatrixPlus__extract_matrix", (DL_FUNC) &BEDMatrixPlus__extract_matrix, 3},
    {"BEDMatrixPlus__extract_columns", (DL_FUNC) &BEDMatrixPlus__extract_columns, 2},
    {"BEDMatrixPlus__new", (DL_FUNC) &BEDMatrixPlus__new, 3},
    {"BEDMatrixPlus__multiply_residuals", (DL_FUNC) &BEDMatrixPlus__multiply_residuals, 5},
    {NULL, NULL, 0}
};

void R_init_BEDMatrixPlus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}