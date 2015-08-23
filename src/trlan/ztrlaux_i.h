#ifndef __ZTRLAN_AUX_H
#define __ZTRLAN_AUX_H
#include "ztrlan.h"
//
////
void trl_print_complex_(trl_info * info, char *title, int size_array,
			trl_dcomplex * array, int inc);
//
// Purpose
// =======
// Print a double precision array for debugging.
//
// Arguments:
// ==========
// info        (input) pointer to the structure trl_info_
//              On entry, points to the data structure to store the information
//              about the eigenvalue problem and the progress of TRLAN
//
// title       (input) character string
//              On entry, specifies the title of the information to be printed.
//
// size_array  (input) integer
//              On entry specifies, the number of doubles to be printed.
//
// array       (input) complex array of length ((size_array-1)*inc+1)
//              On entry, contains the complex values to be printed.
//
// inc         (input) integer
//              On entry, specifies how the index to array should be incremented.
//
//
int ztrl_read_checkpoint(char *filename, int nrow, trl_dcomplex * evec,
			  int lde, int mev, int *j1, trl_dcomplex * base,
			  int ldb, int nbas, int *j2, int nalpha,
			  double *alpha, int nbeta, double *beta);
//
// Purpose
// =======
// Read check-point file
//
// Arguments:
// ==========
// filename   (input) character string
//             On entry, specifies the name of checkpoint file.
//
// nrow       (input) integer
//             On entry, specifies the problem size.
//
// evec       (output) double precision array of dimensioni (lde,j1)
//             On exit, contains the first part of basis vectors stored in checkpoint.
//
// lde         (input) integer
//             On entry, specifies the leading dimension of evec.
//
// mev        (input) integer
//             On entry, specifies the number of eigenvalues converged.
//
// j1         (input) integer
//             On entry, specifies the last column index of evec, that contains a base
//             vector.
//
// base       (output) double precision array of dimension (ldb,nbas)
//             On exit, contains the second part of basis stored in the checkpoint.
//
// ldb        (input) integer
//             On entry, specifies the leading dimension of base.
//
// nbas       (input) integer
//             On entry, specifies the number of columns in base.
//
// j2         (input) integer
//             On entry, specifies the last column index of base, that contains a base
//             vector.
//
// nalpha     (input) integer
//             On entry, specifies the size of alpha
//
// alpha      (output) double precision array of length (nalpha)
//             On exit, contains the alpha values stored in checkpoint.
//
// nbeta      (input) integer
//             On entry, specifies the size of beta.
//
// beta       (output) double precision array of length (nbeta)
//             On exit, contains the beta values stored in checkpoint.
//
int ztrl_write_checkpoint(char *filename, int nrow, double *alpha,
			   double *beta, trl_dcomplex * evec, int lde,
			   int me, trl_dcomplex * base, int ldb, int nb);
//
// Purpose
// =======
// Write a check-point file.
//
// Arguments
// =========
// filename   (input) character string
//             On entry, specifies the name of checkpoint file.
//
// nrow       (input) integer
//             On entry, specifies the problem size.
//
// alpha      (input) double precision array of length (me+nb-1)
//             On entry, contains the alpha values computed so far.
//
// beta       (input) double precision array of length (me+ne-1)
//             On entry, contains the beta values computed so far.
//
// evec       (input) double precision array of dimensioni (lde,me)
//             On entry, contains the first part of basis vectors.
//
// lde        (input) integer
//             On entry, specifies the leading dimension of evec.
//
// me         (input) integer
//             On entry, specifies the last column index of evec, that contains a base
//             vector.
//
// base       (input) double precision array of dimension (ldb,nb)
//             On entry, contains the second part of basis.
//
// ldb        (input) integer
//             On entry, specifies the leading dimension of base.
//
// nb         (input) integer
//             On entry, specifies the last column index of base, that contains a base
//             vector.
//
////
void
ztrl_check_recurrence(ztrl_matprod op, trl_info * info,
                      int nrow, int ncol, trl_dcomplex * v1, int ld1,
                      int m1, trl_dcomplex * v2, int ld2, int m2,
                      int kept, double *alpha, double *beta,
                      trl_dcomplex * wrk, int lwrk, void *lparam);
//
// Purpose
// =======
// Check Lanczos recurrence relation for debug purpose.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information about the
//           eigenvalue problem and the progress of TRLAN.
//
// nrow     (input) integer
//           On entry, specifies the problem size, i.e., the number of ros in v1 and v2.
//
// ncol     (input) integer
//           On entry, specifies the maximum number of eigenvalues that can be stored.
//
// v1       (input) double precision array of dimension (ld1,m1)
//           On entry, contains the first part of basis.
//
// ld1      (input) integer
//           On entry, specifies the leading dimension of v1.
//
// m1       (input) integer
//           On entry, specifies the last column index of v1 that contains a base vector.
//
// v2       (input) double precision array of dimension (ld2,m2)
//           On entry, contains the second part of basis.
//
// ld2      (input) integer
//           On entry, specifies the leading dimension of v2.
//
// m2       (input) integer
//           On entry, specifies the last column index of v2 that contains a base vector.
//
// kept     (input) integer
//           On entry, specifies the number of basis kept at the last restart.
//
// alpha    (input) integer
//           On entry, contains the alpha values computed so far.
//
// beta     (input) integer
//           On entry, contains the beta values computed so far.
//
// wrk      (workspace) double precision vector of length (lwrk)
//
// lwrk     (input) integer
//           On entry, specifies the size of workspace.
// 
////
void ztrl_check_orth(trl_info * info, int nrow, trl_dcomplex * v1,
                     int ld1, int j1, trl_dcomplex * v2, int ld2, int j2,
                     trl_dcomplex * wrk, int lwrk, void *lparam);
//
// Purpose:
// ========
// Check orthogonality of the basis.
//
// Arguments:
// ==========
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information about the 
//           eigenvalue problem and the progress of TRLAN.
//
// nrow     (input) integer
//           On entry, specifies the problem size, i.e., the number of rows in v1 and v2.
//
// v1       (input) double precision array of diimension (ld1,j1)
//           On entry, contains the first part of the basis.
//
// ld1      (input) integer
//           On entry, specifies the leading diimension of v1.
//
// j1       (input) integer
//           On entry, specifies the last column index of v1, containing the basis.
//
// v2       (input) double precision array of dimension (ld2,j2)
//           On entry, contains the second part of the basis.
//
// ld2      (input) integer
//           On entry, specifies the leading dimension of v2.
//
// j2       (input) integer
//           On entry, specifies the last column index of v2, containing the basis.
//
// wrk      (workspace) double precision array of length (lwrk)
//
// lwrk     (input) integer
//           On entry, specifies the size of workspace.
//

#endif
