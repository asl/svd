#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ematmul(SEXP, SEXP, SEXP);
extern SEXP ematmul_unchecked(SEXP, SEXP, SEXP);
extern SEXP extmat_ncol(SEXP);
extern SEXP extmat_nrow(SEXP);
extern SEXP initialize_rextmat(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP is_extmat(SEXP);
extern SEXP propack_svd(SEXP, SEXP, SEXP);
extern SEXP trlan_eigen(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP trlan_svd(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ztrlan_eigen(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ztrlan_svd(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ematmul_",           (DL_FUNC) &ematmul,            3},
    {"ematmul_unchecked",  (DL_FUNC) &ematmul_unchecked,  3},
    {"extmat_ncol",        (DL_FUNC) &extmat_ncol,        1},
    {"extmat_nrow",        (DL_FUNC) &extmat_nrow,        1},
    {"initialize_rextmat", (DL_FUNC) &initialize_rextmat, 5},
    {"is_extmat",          (DL_FUNC) &is_extmat,          1},
    {"propack_svd",        (DL_FUNC) &propack_svd,        3},
    {"trlan_eigen",        (DL_FUNC) &trlan_eigen,        5},
    {"trlan_svd",          (DL_FUNC) &trlan_svd,          5},
    {"ztrlan_eigen",       (DL_FUNC) &ztrlan_eigen,       5},
    {"ztrlan_svd",         (DL_FUNC) &ztrlan_svd,         5},
    {NULL, NULL, 0}
};

void R_init_svd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
