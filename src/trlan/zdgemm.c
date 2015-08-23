/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <stdio.h>

#include "trl_map.h"
#include "trlan.h"

/* Subroutine */ int zdgemm_(int m, int n, int k, trl_dcomplex * a,
			     int lda,
			     double *b, int ldb, trl_dcomplex * c, int ldc)
{
/*  Purpose
    =======

    ZDGEMM  performs the matrix-matrix operations

       C := A*B

    Parameters
    ==========

    M      - INTEGER.
             On entry,  M  specifies  the number  of rows  of the  matrix
             op( A )  and of the  matrix  C.  M  must  be at least  zero.
             Unchanged on exit.

    N      - INTEGER.
             On entry,  N  specifies the number  of columns of the matrix
             op( B ) and the number of columns of the matrix C. N must be
             at least zero.
             Unchanged on exit.

    K      - INTEGER.
             On entry,  K  specifies  the number of columns of the matrix
             op( A ) and the number of rows of the matrix op( B ). K must
             be at least  zero.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
             part of the array  A  must contain the matrix  A,  otherwise
             the leading  k by m  part of the array  A  must contain  the
             matrix A.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then
             LDA must be at least  max( 1, m ), otherwise  LDA must be at
             least  max( 1, k ).
             Unchanged on exit.

    B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
             part of the array  B  must contain the matrix  B,  otherwise
             the leading  n by k  part of the array  B  must contain  the
             matrix B.
             Unchanged on exit.

    LDB    - INTEGER.
             On entry, LDB specifies the first dimension of B as declared
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then
             LDB must be at least  max( 1, k ), otherwise  LDB must be at
             least  max( 1, n ).
             Unchanged on exit.

    C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
             Before entry, the leading  m by n  part of the array  C must
             contain the matrix  C,  except when  beta  is zero, in which
             case C need not be set on entry.
             On exit, the array  C  is overwritten by the  m by n  matrix

             ( alpha*op( A )*op( B ) + beta*C ).

    LDC    - INTEGER.
             On entry, LDC specifies the first dimension of C as declared
             in  the  calling  (sub)  program.   LDC  must  be  at  least
             max( 1, m ).
             Unchanged on exit.   */

    int i1, i2, i3;
    double register tr, ti;

    for (i1 = 0; i1 < n; i1++) {
	for (i2 = 0; i2 < m; i2++) {
	    tr = 0.0;
	    ti = 0.0;
	    for (i3 = 0; i3 < k; i3++) {
		tr += a[i3 * lda + i2].r * b[i1 * ldb + i3];
		ti += a[i3 * lda + i2].i * b[i1 * ldb + i3];
	    }
	    c[i1 * ldc + i2].r = tr;
	    c[i1 * ldc + i2].i = ti;
	}
    }
    return 0;
}
/* Subroutine */ int zdgemm2_(int m, int n, int k, trl_dcomplex * a,
			      int lda,
			      double *b, int ldb, trl_dcomplex * c,
			      int ldc)
{
/*  Purpose
    =======

    ZDGEMM  performs the matrix-matrix operations

       C := A*B+C

    Parameters
    ==========

    M      - INTEGER.
             On entry,  M  specifies  the number  of rows  of the  matrix
             op( A )  and of the  matrix  C.  M  must  be at least  zero.
             Unchanged on exit.

    N      - INTEGER.
             On entry,  N  specifies the number  of columns of the matrix
             op( B ) and the number of columns of the matrix C. N must be
             at least zero.
             Unchanged on exit.

    K      - INTEGER.
             On entry,  K  specifies  the number of columns of the matrix
             op( A ) and the number of rows of the matrix op( B ). K must
             be at least  zero.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
             part of the array  A  must contain the matrix  A,  otherwise
             the leading  k by m  part of the array  A  must contain  the
             matrix A.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then
             LDA must be at least  max( 1, m ), otherwise  LDA must be at
             least  max( 1, k ).
             Unchanged on exit.

    B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
             part of the array  B  must contain the matrix  B,  otherwise
             the leading  n by k  part of the array  B  must contain  the
             matrix B.
             Unchanged on exit.

    LDB    - INTEGER.
             On entry, LDB specifies the first dimension of B as declared
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then
             LDB must be at least  max( 1, k ), otherwise  LDB must be at
             least  max( 1, n ).
             Unchanged on exit.

    C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
             Before entry, the leading  m by n  part of the array  C must
             contain the matrix  C,  except when  beta  is zero, in which
             case C need not be set on entry.
             On exit, the array  C  is overwritten by the  m by n  matrix

             ( alpha*op( A )*op( B ) + beta*C ).

    LDC    - INTEGER.
             On entry, LDC specifies the first dimension of C as declared
             in  the  calling  (sub)  program.   LDC  must  be  at  least
             max( 1, m ).
             Unchanged on exit.   */

    int i1, i2, i3;
    double register tr, ti;

    for (i1 = 0; i1 < n; i1++) {
	for (i2 = 0; i2 < m; i2++) {
	    tr = c[i1 * ldc + i2].r;
	    ti = c[i1 * ldc + i2].i;
	    for (i3 = 0; i3 < k; i3++) {
		tr += a[i3 * lda + i2].r * b[i1 * ldb + i3];
		ti += a[i3 * lda + i2].i * b[i1 * ldb + i3];
	    }
	    c[i1 * ldc + i2].r = tr;
	    c[i1 * ldc + i2].i = ti;
	}
    }
    return 0;
}
