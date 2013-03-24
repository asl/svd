#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>

#include "extmat.h"


typedef struct {
  SEXP fcall;
  SEXP tfcall;
  SEXP rho;
  unsigned n;
  unsigned m;
} rext_matrix;

static void extmat_matmul(double* out,
                          const double* v,
                          const void* matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  SEXP rho, rV, res;
  unsigned n, m;
  PROTECT_INDEX ipx;

  /* Grab the matrix dimensions */
  n = e->n;
  m = e->m;

  /* Grab the environment we're going to evaluate function in */
  rho = e->rho;

  /* Allocate the memory to call R code and prepare the input*/
  PROTECT(rV = allocVector(REALSXP, m));
  Memcpy(REAL(rV), v, m);

  /* Call the actual function */
  SETCADR(e->fcall, rV);
  PROTECT_WITH_INDEX(res = eval(e->fcall, rho), &ipx);
  REPROTECT(res = coerceVector(res, REALSXP), ipx);

  /* Prepare the output */
  Memcpy(out, REAL(res), n);

  UNPROTECT(2);
}

static void extmat_tmatmul(double* out,
                           const double* v,
                           const void* matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  SEXP rho, rV, res;
  unsigned n, m;
  PROTECT_INDEX ipx;

  /* Grab the matrix dimensions */
  n = e->n;
  m = e->m;

  /* Grab the environment we're going to evaluate function in */
  rho = e->rho;

  /* Allocate the memory to call R code and prepare the input*/
  PROTECT(rV = allocVector(REALSXP, n));
  Memcpy(REAL(rV), v, n);

  /* Call the actual function */
  SETCADR(e->tfcall, rV);
  PROTECT_WITH_INDEX(res = eval(e->tfcall, rho), &ipx);
  REPROTECT(res = coerceVector(res, REALSXP), ipx);

  /* Prepare the output */
  Memcpy(out, REAL(res), m);

  UNPROTECT(2);
}

static unsigned extmat_nrow(const void *matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  return e->n;
}

static unsigned extmat_ncol(const void *matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  return e->m;
}

static void extmat_finalizer(SEXP ptr) {
  ext_matrix *e;
  rext_matrix *re;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  e = R_ExternalPtrAddr(ptr);
  if (!e)
    return;

  if (strcmp(e->type, "external matrix from R"))
    return;

  re = (rext_matrix*)e->matrix;

  R_ReleaseObject(re->fcall);
  R_ReleaseObject(re->tfcall);
  R_ReleaseObject(re->rho);

  Free(re);
  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_extmat(SEXP f, SEXP tf, SEXP n, SEXP m, SEXP rho) {
  ext_matrix *e;
  rext_matrix *re;
  SEXP emat;

  /* Allocate memory */
  re = Calloc(1, rext_matrix);

  R_PreserveObject(re->fcall = lang2(f, R_NilValue));
  R_PreserveObject(re->tfcall = lang2(tf, R_NilValue));
  re->n = asInteger(n);
  re->m = asInteger(m);
  R_PreserveObject(re->rho = rho);

  /* Create external matrix envelope */
  e = Calloc(1, ext_matrix);
  e->type = "external matrix from R";
  e->mulfn = extmat_matmul;
  e->tmulfn = extmat_tmatmul;
  e->ncol = extmat_ncol;
  e->nrow = extmat_nrow;

  e->matrix = re;

  /* Make an external pointer envelope */
  emat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue);
  R_RegisterCFinalizer(emat, extmat_finalizer);

  return emat;
}

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
