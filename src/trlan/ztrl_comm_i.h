
#ifndef __ZTRLAN_COMM_H
#define __ZTRLAN_COMM_H

//
void ztrl_g_dot_(int mpicom, int nrow, trl_dcomplex * v1, int ld1, int m1,
		 trl_dcomplex * v2, int ld2, int m2, trl_dcomplex * rr,
		 trl_dcomplex * wrk);
//
// Purpose:
// ========
// Implements a distributed version of BLAS routine dgemv, which is used to compute
// dot-products by TRLAN, i.e., wrk = [V1, V2]'*rr.
//
// Arguments:
// ==========
// mpicom     (input) integer
//             On entry, specifies MPI communicator.
//
// nrow       (input) integer
//             On entry, specifies, the number of rows on the local processor.
//
// v1         (input) double precision array of dimension (ld1,m1)
//             On entry, contains the first part of the matrix.
//
// ld1        (input) integer
//             On entry, specifies the leading dimension of v1.
//
// m1         (input) integer
//             On entry, specifies the number of columns in v1.
//
// v2         (input) double precision array of dimension (ld2,m2)
//             On entry, contains the second part of the matrix.
//
// ld2        (input) integer
//             On entry, specifies the leading dimension of v2.
//
// m2         (input) integer
//             On entry, specifies the number of columns in v2.
//
// rr         (input) double precision array of length (nrow)
//             On entry, contains the vector to be multiplied.
//
// wrk        (output) double precision array of length (m1+m2)
//             On exit, contains th results of this operation.  !! size not checked !!
//
////
void ztrl_g_sum(int mpicom, int nelm, trl_dcomplex * x, trl_dcomplex * y);
//
// Purpose:
// ========
// Performs global sum in the parallel environment, nothing is done here.
//
// Arguments:
// ==========
// mpicom    (input) integer
//            On entry, specifies the MPI communicator.
//
// nelm      (input) integer
//            On entry, specifies the number of elements in x and y.
//
// x         (input/output) double precision array of dimension nelm
//            On entry, contains the resulting vector on this processor.
//            On exit, contain the resulting vector of global sum.
//
// y         (workspace) double precision array of dimensioni nelm
//

#endif
 
