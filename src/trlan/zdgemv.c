/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <stdio.h>

#include "ztrlan.h"
#include "trlan.h"

/* Subroutine */ int zdgemv_(int m, int n, double alpha, trl_dcomplex * a,
			     int lda,
			     double *x, double beta, trl_dcomplex * y)
{


/*  Purpose
    =======

    ZGEMV  performs one of the matrix-vector operations

       y := beta * y - alpha * A * x

    where alpha and beta are scalars, x and y are vectors and A is an
    m by n matrix.

    Parameters
    ==========


    M      - INTEGER.
             On entry, M specifies the number of rows of the matrix A.
             M must be at least zero.
             Unchanged on exit.

    N      - INTEGER.
             On entry, N specifies the number of columns of the matrix A.
             N must be at least zero.
             Unchanged on exit.

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
             Before entry, the leading m by n part of the array A must
             contain the matrix of coefficients.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared

             in the calling (sub) program. LDA must be at least
             max( 1, m ).
             Unchanged on exit.

    X      - COMPLEX*16       array of DIMENSION at least
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
             and at least
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
             Before entry, the incremented array X must contain the
             vector x.
             Unchanged on exit.

    Y      - COMPLEX*16       array of DIMENSION at least
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
             and at least
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
             Before entry with BETA non-zero, the incremented array Y
             must contain the vector y. On exit, Y is overwritten by the

             updated vector y.


    Level 2 Blas routine.

    -- Written on 22-October-1986.
       Jack Dongarra, Argonne National Lab.
       Jeremy Du Croz, Nag Central Office.
       Sven Hammarling, Nag Central Office.
       Richard Hanson, Sandia National Labs.



       Test the input parameters.


   Parameter adjustments
       Function Body */
    int i1, i2;
    double register tr, ti;

    for (i1 = 0; i1 < m; i1++) {
	tr = 0.0;
	ti = 0.0;
	for (i2 = 0; i2 < n; i2++) {
	    tr += a[i2 * lda + i1].r * x[i2];
	    ti += a[i2 * lda + i1].i * x[i2];
	}
	y[i1].r = beta * y[i1].r + alpha * tr;
	y[i1].i = beta * y[i1].i + alpha * ti;
    }
    return 0;
}
