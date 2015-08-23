#ifndef __ZTRLAN_I_H__
#define __ZTRLAN_I_H__

#include "ztrlan.h"

void trl_zaxpy(int n, trl_dcomplex za, trl_dcomplex * zx, int incx,
               trl_dcomplex * zy, int incy);
void trl_zgemv(char *trans, int m, int n, trl_dcomplex alpha,
               trl_dcomplex * a, int lda, trl_dcomplex * x, int incx,
               trl_dcomplex beta, trl_dcomplex * y, int incy);
void trl_zgemm(char *transa, char *transb, int m, int n, int k,
               trl_dcomplex alpha, trl_dcomplex * a, int lda,
               trl_dcomplex * b, int ldb, trl_dcomplex beta,
               trl_dcomplex * c, int ldc);
void trl_zscal(int n, trl_dcomplex za, trl_dcomplex * zx, int incx);
void trl_zdscal(int n, double da, trl_dcomplex * zx, int incx);
void trl_zcopy(int n, trl_dcomplex * zx, int incx, trl_dcomplex * zy,
               int incy);
void trl_zdotc(trl_dcomplex * ret_val, int n, trl_dcomplex * zx, int incx,
               trl_dcomplex * zy, int incy);

#endif
