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

static void rextmat_matmul(double* out,
                           const double* v,
                           const void* matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  SEXP rho, rV, res, fcall;
  unsigned n, m;
  PROTECT_INDEX ipx;

  /* Grab the matrix dimensions */
  n = e->n;
  m = e->m;

  /* Grab the environment we're going to evaluate function in */
  rho = R_WeakRefValue(e->rho);

  /* Grab the function */
  fcall = R_WeakRefValue(e->fcall);

  /* Allocate the memory to call R code and prepare the input*/
  PROTECT(rV = allocVector(REALSXP, m));
  Memcpy(REAL(rV), v, m);

  /* Call the actual function */
  SETCADR(fcall, rV);
  PROTECT_WITH_INDEX(res = eval(fcall, rho), &ipx);
  REPROTECT(res = coerceVector(res, REALSXP), ipx);

  /* Prepare the output */
  Memcpy(out, REAL(res), n);

  UNPROTECT(2);
}

static void rextmat_tmatmul(double* out,
                            const double* v,
                            const void* matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  SEXP rho, rV, res, tfcall;
  unsigned n, m;
  PROTECT_INDEX ipx;

  /* Grab the matrix dimensions */
  n = e->n;
  m = e->m;

  /* Grab the environment we're going to evaluate function in */
  rho = R_WeakRefValue(e->rho);

  /* Grab the function */
  tfcall = R_WeakRefValue(e->tfcall);

  /* Allocate the memory to call R code and prepare the input*/
  PROTECT(rV = allocVector(REALSXP, n));
  Memcpy(REAL(rV), v, n);

  /* Call the actual function */
  SETCADR(tfcall, rV);
  PROTECT_WITH_INDEX(res = eval(tfcall, rho), &ipx);
  REPROTECT(res = coerceVector(res, REALSXP), ipx);

  /* Prepare the output */
  Memcpy(out, REAL(res), m);

  UNPROTECT(2);
}

static unsigned rextmat_nrow(const void *matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  return e->n;
}

static unsigned rextmat_ncol(const void *matrix) {
  rext_matrix *e = (rext_matrix*)matrix;

  return e->m;
}

static void rextmat_finalizer(SEXP ptr) {
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

  Free(re);
  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_rextmat(SEXP f, SEXP tf, SEXP n, SEXP m, SEXP rho) {
  ext_matrix *e;
  rext_matrix *re;
  SEXP emat;

  /* Allocate memory */
  re = Calloc(1, rext_matrix);

  re->n = asInteger(n);
  re->m = asInteger(m);

  /* Create external matrix envelope */
  e = Calloc(1, ext_matrix);
  e->type = "external matrix from R";
  e->mulfn = rextmat_matmul;
  e->tmulfn = rextmat_tmatmul;
  e->ncol = rextmat_ncol;
  e->nrow = rextmat_nrow;

  e->matrix = re;

  /* Make an external pointer envelope */
  PROTECT(emat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue));

  /* Attach the fields */
  PROTECT(re->fcall = R_MakeWeakRef(emat, lang2(f, R_NilValue), R_NilValue, 1));
  PROTECT(re->tfcall = R_MakeWeakRef(emat, lang2(tf, R_NilValue), R_NilValue, 1));
  PROTECT(re->rho = R_MakeWeakRef(emat, rho, R_NilValue, 1));

  R_RegisterCFinalizer(emat, rextmat_finalizer);

  UNPROTECT(4);

  return emat;
}

SEXP ematmul_unchecked(SEXP emat, SEXP v, SEXP transposed) {
  SEXP Y = NILSXP;
  R_len_t K, L;
  ext_matrix *e;
  void *matrix;

  /* Grab needed data */
  e = R_ExternalPtrAddr(emat);
  matrix = e->matrix;

  L = (LOGICAL(transposed)[0] ? e->ncol(matrix) : e->nrow(matrix));
  K = (LOGICAL(transposed)[0] ? e->nrow(matrix) : e->ncol(matrix));

  /* Check agains absurd values of inputs */
  if (K != length(v))
    error("invalid length of input vector 'v'");

  /* Allocate output buffer */
  PROTECT(Y = allocVector(REALSXP, L));

  /* Calculate the product */
  if (LOGICAL(transposed)[0])
    e->tmulfn(REAL(Y), REAL(v), matrix);
  else
    e->mulfn(REAL(Y), REAL(v), matrix);

  UNPROTECT(1);

  return Y;
}

SEXP ematmul(SEXP emat, SEXP v, SEXP transposed) {
  SEXP Y = NILSXP, tchk;

  /* Perform a type checking */
  PROTECT(tchk = is_extmat(emat));

  if (LOGICAL(tchk)[0]) {
    Y = ematmul_unchecked(emat, v, transposed);
  } else
    error("pointer provided is not an external matrix");

  UNPROTECT(1);

  return Y;
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

SEXP extmat_nrow(SEXP ptr) {
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

SEXP extmat_ncol(SEXP ptr) {
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