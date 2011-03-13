#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "extmat.h"

SEXP is_extmat(SEXP ptr) {
  SEXP ans;
  ext_matrix *e = NULL;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* object is an external pointer */
  if (TYPEOF(ptr) != EXTPTRSXP)
    LOGICAL(ans)[0] = 0;

  /* tag should be 'external matrix' */
  if (LOGICAL(ans)[0] &&
      R_ExternalPtrTag(ptr) != install("external matrix"))
    LOGICAL(ans)[0] = 0;

  /* pointer itself should not be null */
  if (LOGICAL(ans)[0]) {
    e = R_ExternalPtrAddr(ptr);
    if (!e)
      LOGICAL(ans)[0] = 0;
  }

  /* finally, type should be nonnull */
  if (LOGICAL(ans)[0] && e && e->type == NULL)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(1);

  return ans;
}

SEXP extmat_rows(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_extmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = e->nrow(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not an external matrix");

  UNPROTECT(1);

  return ans;
}

SEXP extmat_cols(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_extmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = e->ncol(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not an external matrix");

  UNPROTECT(1);

  return ans;
}

void R_init_svd(DllInfo *info) {
  (void)info;

  R_RegisterCCallable("svd", "is_extmat", (DL_FUNC)is_extmat);
  R_RegisterCCallable("svd", "extmat_rows", (DL_FUNC)extmat_rows);
  R_RegisterCCallable("svd", "extmat_cols", (DL_FUNC)extmat_cols);
}

