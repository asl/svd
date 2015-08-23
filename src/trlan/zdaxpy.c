/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/
#include <stdio.h>

#include "ztrlan.h"
#include "trlan.h"

/* Subroutine */ int zdaxpy_(int n, double a, trl_dcomplex * zx,
			     trl_dcomplex * zy)
{
/*
    Purpose
    =======
    ZDAXPY performs the axpy operation
       zy = zy + a * zx
    where zy and zx are complex and a is real.

    Parameters
    ==========
    N      - INTEGER.
             On entry,  N  specifies the size of the vectors zx and zy.
             Unchanged on exit.

    A      - DOUBLE
             On entry,  A  specifies the scalar part of the axpy operation.

    ZX     - DOUBLE COMPLEX
             On entry,   ZX  contains the vector to update ZY with in the axpy operation.

    ZY     - DOUBLE COMPLEX
             On entry,   ZY  contains the vector to be updated in the axpy operation.
             On exit,    ZY  contains the update vector.
*/
    int i;
    for (i = 0; i < n; i++) {
	zy[i].r += a * zx[i].r;
	zy[i].i += a * zx[i].i;
    }

    return 0;
}
