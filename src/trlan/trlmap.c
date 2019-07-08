/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R_ext/BLAS.h>
#include <R_ext/Complex.h>

typedef Rcomplex trl_dcomplex;

double trl_ddot(int n, const double *dx, int incx,
                const double *dy, int incy) {
  return F77_CALL(ddot)(&n, dx, &incx, dy, &incy);
}

void trl_dcopy(int n, double *dx, int incx, double *dy, int incy) {

  F77_CALL(dcopy)(&n, dx, &incx, dy, &incy);
}

void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
               double *x, int incx, double beta, double *y, int incy) {
  F77_CALL(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}

void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
               int incy) {
  F77_CALL(daxpy)(&n, &da, dx, &incx, dy, &incy);
}

void trl_dgemm(char *transa, char *transb, int m, int n, int k,
               double alpha, double *a, int lda, double *b, int ldb,
               double beta, double *c, int ldc) {
  F77_CALL(dgemm)(transa, transb,
                  &m, &n, &k,
                  &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

void trl_dscal(int n, double da, double *dx, int incx) {
  F77_CALL(dscal)(&n, &da, dx, &incx);
}

void trl_zaxpy(int n, trl_dcomplex za, trl_dcomplex * zx, int incx,
               trl_dcomplex * zy, int incy) {
  F77_CALL(zaxpy)(&n, &za, zx, &incx, zy, &incy);
}

void trl_zgemv(char *trans, int m, int n, trl_dcomplex alpha,
               trl_dcomplex * a, int lda, trl_dcomplex * x, int incx,
               trl_dcomplex beta, trl_dcomplex * y, int incy) {
  F77_CALL(zgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}

void trl_zgemm(char *transa, char *transb, int m, int n, int k,
	       trl_dcomplex alpha, trl_dcomplex * a, int lda,
	       trl_dcomplex * b, int ldb, trl_dcomplex beta,
	       trl_dcomplex * c, int ldc) {
  F77_CALL(zgemm)(transa, transb,
                  &m, &n, &k,
                  &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

void trl_zscal(int n, trl_dcomplex za, trl_dcomplex * zx, int incx) {
  F77_CALL(zscal)(&n, &za, zx, &incx);
}

void trl_zdscal(int n, double da, trl_dcomplex * zx, int incx) {
  F77_CALL(zdscal)(&n, &da, zx, &incx);
}

void trl_zcopy(int n, trl_dcomplex * zx, int incx, trl_dcomplex * zy,
               int incy) {
  F77_CALL(zcopy)(&n, zx, &incx, zy, &incy);
}

void trl_zdotc(trl_dcomplex * ret_val, int n, trl_dcomplex * zx, int incx,
               trl_dcomplex * zy, int incy) {
  *ret_val = F77_CALL(zdotc)(&n, zx, &incx, zy, &incy);
}
