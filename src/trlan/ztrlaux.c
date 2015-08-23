/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "ztrlcore_i.h"
#include "ztrl_comm_i.h"
#include "ztrlaux_i.h"
#include "ztrlan_i.h"

////
void trl_print_complex_(trl_info * info, char *title, int size_array,
			trl_dcomplex * array, int inc)
{
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
// ..
// .. local scalars ..
    int i;
//
// ..
// .. executable statements ..
    fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
    if (size_array > 1) {
	fprintf(info->log_fp, "\n");
    }
    for (i = 0; i < size_array; i += inc) {
	fprintf(info->log_fp, " %10.7e+%10.7ei", array[i].r, array[i].i);
	if ((i % 4) == 3)
	    fprintf(info->log_fp, "\n");
    }
    if (((size_array - 1) % 4) != 3)
	fprintf(info->log_fp, "\n");
//
// .. end of trl_print_complex_ ..
}

////
int ztrl_write_checkpoint(char *filename, int nrow, double *alpha,
			   double *beta, trl_dcomplex * evec, int lde,
			   int me, trl_dcomplex * base, int ldb, int nb)
{
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
// ..
// .. local variables ..
    int jnd, i, j;
    FILE *io_fp;
//
// ..
// .. executable statements ..
    jnd = me + nb - 1;
    io_fp = fopen(filename, "w");
    if (io_fp == NULL) {
      Rprintf("TRL_WRITE_CHECKPOINT: failed to open file: %s.\n",
	       filename);
	return -221;
    }
    if (fwrite(&nrow, sizeof(nrow), 1, io_fp) < 1) {
	return close_file(io_fp, -223, -222);
    }
    if (fwrite(&jnd, sizeof(jnd), 1, io_fp) < 1) {
	return close_file(io_fp, -223, -222);
    }

    for (i = 0; i < jnd; i++) {
	if (fwrite(&alpha[i], sizeof(alpha[i]), 1, io_fp) < 1) {
	    return close_file(io_fp, -223, -222);
	}
    }
    for (i = 0; i < jnd; i++) {
	if (fwrite(&beta[i], sizeof(beta[i]), 1, io_fp) < 1) {
	    return close_file(io_fp, -223, -222);
	}
    }
    for (i = 0; i < me; i++) {
	for (j = 0; j < nrow; j++) {
	    if (fwrite
		(&evec[i * lde + j], sizeof(evec[i * lde + j]), 1,
		 io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	    }
	}
    }
    for (i = 0; i < nb; i++) {
	for (j = 0; j < nrow; j++) {
	    if (fwrite
		(&base[i * ldb + j], sizeof(base[i * ldb + j]), 1,
		 io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	    }
	}
    }
    return close_file(io_fp, 0, -223);
//
// .. end of trl_write_checkpoint ..
//
}

////
int ztrl_read_checkpoint(char *filename, int nrow, trl_dcomplex * evec,
			  int lde, int mev, int *j1, trl_dcomplex * base,
			  int ldb, int nbas, int *j2, int nalpha,
			  double *alpha, int nbeta, double *beta)
{
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
// ..
// .. local variables ..
    int i, j;
    FILE *io_fp;
//
// ..
// .. executable statements ..
    if (lde < nrow || ldb < nrow) {
      Rprintf("TRL_READ_CHECKPOINT: leading dimensions too small.\n");
	return -211;
    }
    // open file
    io_fp = fopen(filename, "r");
    if (io_fp == NULL) {
      Rprintf
	    ("TRL_READ_CHECKPOINT: failed to open check-point file %s.\n",
	     filename);
	return -212;
    }
    // read size information
    if (fread(j1, sizeof(*j1), 1, io_fp) <= 0) {
	return close_file(io_fp, -215, -216);
    }
    if (fread(j2, sizeof(*j2), 1, io_fp) <= 0) {
	return close_file(io_fp, -215, -216);
    }
    if (*j1 != nrow) {
      Rprintf("TRL_READ_CHECKPOINT: Nrow mismatch.\n");
	return -213;
    }
    if (*j2 > imin2(nalpha, imin2(nbeta, mev + nbas - 1))) {
      Rprintf("TRL_READ_CHECKPOINT: MAXLAN too small.");
	return -214;
    }
    // can continue read all data
    for (i = 0; i < *j2; i++) {
	if (fread(&alpha[i], sizeof(alpha[i]), 1, io_fp) <= 0) {
	    return close_file(io_fp, -215, -216);
	}
    }
    for (i = 0; i < *j2; i++) {
	if (fread(&beta[i], sizeof(beta[i]), 1, io_fp) <= 0) {
	    return close_file(io_fp, -215, -216);
	}
    }
    *j1 = imin2(mev, *j2);
    *j2 = *j2 - *j1;
    if (*j1 < mev) {
	for (i = 0; i <= *j1; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&evec[i * lde + j], sizeof(evec[i * lde + j]), 1,
		     io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
    } else {
	for (i = 0; i < *j1; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&evec[i * lde + j], sizeof(evec[i * lde + j]), 1,
		     io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
	for (i = 0; i < *j2; i++) {
	    for (j = 0; j < nrow; j++) {
		if (fread
		    (&base[i * ldb + j], sizeof(base[i * ldb + j]), 1,
		     io_fp) <= 0) {
		    return close_file(io_fp, -215, -216);
		}
	    }
	}
    }
    return close_file(io_fp, 0, -215);
//
// .. end of trl_read_checkpoint ..
//
}

////
void ztrl_check_orth(trl_info * info, int nrow, trl_dcomplex * v1,
		      int ld1, int j1, trl_dcomplex * v2, int ld2, int j2,
                     trl_dcomplex * wrk, int lwrk, void *lparam)
{
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
//
// ..
// .. local parameters ..
    double one = 1.0, zero = 0.0;
    long c__1 = 1;
//
// ..
// .. local variables
    int i, j, k, jnd;
    double nrmfro, nrminf;
    trl_dcomplex tmp;
//
// ..
// .. executable statements ..
    jnd = j1 + j2;
    nrmfro = zero;
    nrminf = zero;
    if (jnd <= 0)
	return;
    if (lwrk < (jnd + jnd)) {
	fprintf(info->log_fp, "TRL_CHECK_ORTH: workspace too small.\n");
	return;
    }
    fprintf(info->log_fp,
	    "TRL_CHECK_ORTH: check orthogonality of %d basis vectors.\n",
	    jnd);
    //
    // check orthognality of the basis vectors
    //
    for (i = 0; i < j1; i++) {
	ztrl_g_dot_(info->mpicom, nrow, v1, ld1, i + 1, v2, ld2, 0,
		    &v1[i * ld1], wrk);
	wrk[i].r = wrk[i].r - one;
	if (info->verbose > 7) {
	    fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
		    i + 1);
	    for (j = 0; j <= i; j++) {
		fprintf(info->log_fp, " %10.3e + i * %10.3e, ", wrk[j].r,
			wrk[j].i);
		if ((j % 5) == 4)
		    fprintf(info->log_fp, "\n");
	    }
	    if ((i % 5) != 4)
		fprintf(info->log_fp, "\n");
	}
	trl_zdotc(&tmp, i, wrk, c__1, wrk, c__1);
	nrmfro += (2 * tmp.r + wrk[i].r * wrk[i].r + wrk[i].i * wrk[i].i);
	if (i == 0) {
	    wrk[i + 1].r = fabs(wrk[i].r);
	} else {
	    wrk[i + 1].r = fmax2(wrk[i].r, wrk[i - 1].r);
	}
	nrminf = fmax2(nrminf, wrk[i + 1].r);
    }
    for (i = 0; i < j2; i++) {
	j = j1 + i;
	ztrl_g_dot_(info->mpicom, nrow, v1, ld1, j1, v2, ld2, i + 1,
		    &v2[i * ld2], wrk);
	wrk[j].r = wrk[j].r - one;
	if (info->verbose > 7) {
	    fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
		    j + 1);
	    for (k = 0; k <= j; k++) {
		fprintf(info->log_fp, " %10.3e + i * %10.3e, ", wrk[k].r,
			wrk[k].i);
		if ((k % 5) == 4)
		    fprintf(info->log_fp, "\n");
	    }
	    if ((j % 5) != 4)
		fprintf(info->log_fp, "\n");
	}
	trl_zdotc(&tmp, j, wrk, c__1, wrk, c__1);
	nrmfro += (2 * tmp.r + wrk[j].r * wrk[j].r + wrk[j].i * wrk[j].i);
	nrminf = fmax2(nrminf, fabs(wrk[j].r));
    }
    fprintf(info->log_fp,
	    "Frobenius norm of orthogonality level %10i %4i  %14.5e\n",
	    info->matvec, jnd, sqrt(nrmfro));
    fprintf(info->log_fp,
	    "Maximum abs. value of orthogonality level is  %14.5e\n",
	    nrminf);
//
// .. end of trl_check_orth ..
//
}

////
void
ztrl_check_recurrence(ztrl_matprod op, trl_info * info,
                      int nrow, int ncol, trl_dcomplex * v1, int ld1,
                      int m1, trl_dcomplex * v2, int ld2, int m2,
                      int kept, double *alpha, double *beta,
                      trl_dcomplex * wrk, int lwrk, void *lparam)
{
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
//
// ..
// .. local parameters ..
    int i__1 = 1;
    long c__1 = 1;
    double zero = 0.0, one = 1.0;
//
// ..
// .. local variables ..
    int i, ii, j, j1, j2, jnd, mv1, ldqkp1;
    char title[TRLAN_STRING_LEN];
    trl_dcomplex *aq, *qkp1, *cs, *alf, *bet;
//
// ..
// .. executable statements ..
    mv1 = m1;
    if (m2 > 0) {
	j2 = m2 - 1;
	j1 = m1;
    } else {
	j2 = 0;
	j1 = m1 - 1;
    }
    jnd = j1 + j2;
    if (lwrk < jnd * 4 + imax2(jnd * 4, nrow)) {
	fprintf(info->log_fp,
		"TRL_CHECK_RECURRENCE: not enough workspace.\n");
	return;
    }
    if (lwrk >= jnd * 4 + nrow) {
	aq = &wrk[lwrk - nrow];
    } else if (lwrk >= jnd * 4) {
      aq = Calloc(nrow, trl_dcomplex);
	if (aq == NULL) {
	    fprintf(info->log_fp,
		    "TRL_CHECK_RECURRENCE: failed to allcoate workspace.\n");
	    return;
	}
    }
    memset(wrk, 0, 4 * jnd * sizeof(trl_dcomplex));
    cs = &wrk[jnd];
    alf = &wrk[2 * jnd];
    bet = &wrk[3 * jnd];
    //
    // first type of relation
    // A q_i = Alpha_i q_i + Beta_i q_{k+1}
    if (kept < ncol) {
	qkp1 = &v1[kept * ld1];
	ldqkp1 = ld1;
    } else {
	qkp1 = &v2[(kept - j1) * ld2];
	ldqkp1 = ld2;
    }

    // alf(i) = v1(:,i)' * A * v1(:,i)
    // aq = A * v1(:,i) - alpha(i) * v1(:,i)
    // bet(i) = aq' * aq
    // cs(i) = aq' *
    for (i = 0; i < imin2(j1, kept); i++) {
      op(&nrow, &i__1, &v1[i * ld1], &ld1, aq, &nrow, lparam);
	//printf( "1: alpha[%d]\n",i );
	for (ii = 0; ii < nrow; ii++) {
	    alf[i].r +=
		(aq[ii].r * v1[i * ld1 + ii].r +
		 aq[ii].i * v1[i * ld1 + ii].i);
	    alf[i].i +=
		(aq[ii].i * v1[i * ld1 + ii].r -
		 aq[ii].r * v1[i * ld1 + ii].i);

	    aq[ii].r -= (alpha[i] * v1[i * ld1 + ii].r);
	    aq[ii].i -= (alpha[i] * v1[i * ld1 + ii].i);
	    bet[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

	    cs[i].r += (aq[ii].r * qkp1[ii].r + aq[ii].i * qkp1[ii].i);
	    cs[i].i += (aq[ii].r * qkp1[ii].i - aq[ii].i * qkp1[ii].r);

	    aq[ii].r -= (beta[i] * qkp1[ii].r);
	    aq[ii].i -= (beta[i] * qkp1[ii].i);

	    wrk[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	}
    }
    for (i = 0; i < (kept - j1); i++) {
	j = i + j1;
	op(&nrow, &i__1, &v2[i * ld2], &ld2, aq, &nrow, lparam);
	//printf( "2: alpha[%d]\n",j );
	for (ii = 0; ii < nrow; ii++) {
	    alf[j].r +=
		(aq[ii].r * v2[i * ld2 + ii].r +
		 aq[ii].i * v2[i * ld2 + ii].i);
	    alf[j].i +=
		(aq[ii].i * v2[i * ld2 + ii].r -
		 aq[ii].r * v2[i * ld2 + ii].i);

	    aq[ii].r -= (alpha[j] * v2[i * ld2 + ii].r);
	    aq[ii].i -= (alpha[j] * v2[i * ld2 + ii].i);
	    bet[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

	    cs[j].r += (aq[ii].r * qkp1[ii].r + aq[ii].i * qkp1[ii].i);
	    cs[j].i += (aq[ii].r * qkp1[ii].i - aq[ii].i * qkp1[ii].r);

	    aq[ii].r -= (beta[j] * qkp1[ii].r);
	    aq[ii].i -= (beta[j] * qkp1[ii].i);

	    wrk[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	}
    }
    //
    // the (k+1)st base vector need to orthogonalize against all previous
    // vectors
    if (jnd > kept) {
      op(&nrow, &i__1, qkp1, &ldqkp1, aq, &nrow, lparam);
	trl_zdotc(&(alf[kept]), nrow, aq, c__1, qkp1, c__1);
	zdaxpy_(nrow, -alpha[kept], qkp1, aq);
	for (i = 0; i < imin2(j1, kept); i++) {
	    zdaxpy_(nrow, -beta[i], &v1[i * ld1], aq);
	}
	for (i = 0; i < kept - j1; i++) {
	    j = j1 + i;
	    zdaxpy_(nrow, -beta[j], &v2[i * ld2], aq);
	}
	trl_zdotc(&(bet[kept]), nrow, aq, c__1, aq, c__1);
	if (kept + 2 <= j1) {
	    trl_zdotc(&(cs[kept]), nrow, aq, c__1, &v1[(kept + 1) * ld1],
		      c__1);
	    zdaxpy_(nrow, -beta[kept], &v1[(kept + 1) * ld1], aq);
	} else {
	    trl_zdotc(&(cs[kept]), nrow, aq, c__1,
		      &v2[(kept + 1 - j1) * ld2], c__1);
	    zdaxpy_(nrow, -beta[kept], &v2[(kept + 1 - j1) * ld2], aq);
	}
	trl_zdotc(&(wrk[kept]), nrow, aq, c__1, aq, c__1);
    }
    //
    // the third kind of relation -- normal three term recurrence
    // depending the fact that if the lower-bound of loop is less than
    // upper bound, the look should not be executed
    for (i = kept + 1; i < j1; i++) {
      op(&nrow, &i__1, &v1[ld1 * i], &ld1, aq, &nrow, lparam);
	if (i < (mv1 - 1)) {
	    for (ii = 0; ii < nrow; ii++) {
		alf[i].r +=
		    (aq[ii].r * v1[i * ld1 + ii].r +
		     aq[ii].i * v1[i * ld1 + ii].i);
		alf[i].i +=
		    (aq[ii].i * v1[i * ld1 + ii].r -
		     aq[ii].r * v1[i * ld1 + ii].i);

		aq[ii].r -=
		    (alpha[i] * v1[i * ld1 + ii].r +
		     beta[i - 1] * v1[(i - 1) * ld1 + ii].r);
		aq[ii].i -=
		    (alpha[i] * v1[i * ld1 + ii].i +
		     beta[i - 1] * v1[(i - 1) * ld1 + ii].i);
		bet[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

		cs[i].r +=
		    (aq[ii].r * v1[(i + 1) * ld1 + ii].r +
		     aq[ii].i * v1[(i + 1) * ld1 + ii].i);
		cs[i].i +=
		    (aq[ii].r * v1[(i + 1) * ld1 + ii].i -
		     aq[ii].i * v1[(i + 1) * ld1 + ii].r);

		aq[ii].r -= beta[i] * v1[(i + 1) * ld1 + ii].r;
		aq[ii].i -= beta[i] * v1[(i + 1) * ld1 + ii].i;

		wrk[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	    }
	} else {
	    for (ii = 0; ii < nrow; ii++) {
		alf[i].r +=
		    (aq[ii].r * v1[i * ld1 + ii].r +
		     aq[ii].i * v1[i * ld1 + ii].i);
		alf[i].i +=
		    (aq[ii].i * v1[i * ld1 + ii].r -
		     aq[ii].r * v1[i * ld1 + ii].i);

		aq[ii].r -=
		    (alpha[i] * v1[i * ld1 + ii].r +
		     beta[i - 1] * v1[(i - 1) * ld1 + ii].r);
		aq[ii].i -=
		    (alpha[i] * v1[i * ld1 + ii].i +
		     beta[i - 1] * v1[(i - 1) * ld1 + ii].i);
		bet[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

		cs[i].r += (aq[ii].r * v2[ii].r + aq[ii].i * v2[ii].i);
		cs[i].i += (aq[ii].r * v2[ii].i - aq[ii].i * v2[ii].r);

		aq[ii].r -= beta[i] * v2[ii].r;
		aq[ii].i -= beta[i] * v2[ii].i;
		wrk[i].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	    }
	}
    }
    for (i = imax2(0, kept - j1 + 1); i < j2; i++) {
	j = i + j1;
	op(&nrow, &i__1, &v2[i * ld2], &ld2, aq, &nrow, lparam);
	if (i > 0) {
	    for (ii = 0; ii < nrow; ii++) {
		alf[j].r +=
		    (aq[ii].r * v2[i * ld2 + ii].r +
		     aq[ii].i * v2[i * ld2 + ii].i);
		alf[j].i +=
		    (aq[ii].i * v2[i * ld2 + ii].r -
		     aq[ii].r * v2[i * ld2 + ii].i);

		aq[ii].r -=
		    (beta[j - 1] * v2[(i - 1) * ld2 + ii].r +
		     alpha[j] * v2[i * ld2 + ii].r);
		aq[ii].i -=
		    (beta[j - 1] * v2[(i - 1) * ld2 + ii].i +
		     alpha[j] * v2[i * ld2 + ii].i);
		bet[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

		cs[j].r +=
		    (aq[ii].r * v2[(i + 1) * ld2 + ii].r +
		     aq[ii].i * v2[(i + 1) * ld2 + ii].i);
		cs[j].i +=
		    (aq[ii].r * v2[(i + 1) * ld2 + ii].i -
		     aq[ii].i * v2[(i + 1) * ld2 + ii].r);

		aq[ii].r -= beta[j] * v2[(i + 1) * ld2 + ii].r;
		aq[ii].i -= beta[j] * v2[(i + 1) * ld2 + ii].i;
		wrk[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	    }
	} else {
	    for (ii = 0; ii < nrow; ii++) {
		alf[j].r += (aq[ii].r * v2[ii].r + aq[ii].i * v2[ii].i);
		alf[j].i += (aq[ii].i * v2[ii].r - aq[ii].r * v2[ii].i);

		aq[ii].r -=
		    (beta[j - 1] * v1[(j1 - 1) * ld1 + ii].r +
		     alpha[j] * v2[ii].r);
		aq[ii].i -=
		    (beta[j - 1] * v1[(j1 - 1) * ld1 + ii].i +
		     alpha[j] * v2[ii].i);
		bet[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);

		cs[j].r +=
		    (aq[ii].r * v2[ld2 + ii].r +
		     aq[ii].i * v2[ld2 + ii].i);
		cs[j].i +=
		    (aq[ii].r * v2[ld2 + ii].i -
		     aq[ii].i * v2[ld2 + ii].r);

		aq[ii].r -= beta[j] * v2[ld2 + ii].r;
		aq[ii].i -= beta[j] * v2[ld2 + ii].i;
		wrk[j].r += (aq[ii].r * aq[ii].r + aq[ii].i * aq[ii].i);
	    }
	}
    }
    //
    ztrl_g_sum(info->mpicom, jnd * 4, wrk, &wrk[jnd * 4]);
    aq[0].r = zero;
    for (ii = 0; ii < jnd; ii++) {
	aq[0].r += wrk[ii].r;
    }
    aq[0].r = sqrt(aq[0].r);
    for (ii = 0; ii < jnd; ii++) {
	wrk[ii].r = sqrt(wrk[ii].r);
    }
    for (ii = 0; ii < jnd; ii++) {
	if (bet[ii].r > zero) {
	    if (beta[ii] < zero) {
		bet[ii].r = -sqrt(bet[ii].r);
	    } else {
		bet[ii].r = sqrt(bet[ii].r);
	    }
	    cs[ii].r = cs[ii].r / bet[ii].r;
	    cs[ii].i = cs[ii].i / bet[ii].r;
	} else {
	    bet[ii].r = zero;
	}
    }
    strcpy(title, "Alpha computed by TRLAN ..");
    trl_print_real(info, title, jnd, alpha, 1);
    strcpy(title, "Alpha computed explicitly in TRL_CHECK_RECURRENCE ..");
    trl_print_complex_(info, title, jnd, alf, 1);
    strcpy(title, "Differences in alpha ..");
    for (ii = 0; ii < jnd; ii++) {
	alf[ii].r -= alpha[ii];
    }
    trl_print_complex_(info, title, jnd, alf, 1);
    strcpy(title, "Beta computed by TRLAN ..");
    trl_print_real(info, title, jnd, beta, 1);
    strcpy(title, "Beta computed explicitly in TRL_CHECK_RECURRENCE ..");
    trl_print_complex_(info, title, jnd, bet, 1);
    strcpy(title, "Differences in beta ..");
    for (ii = 0; ii < jnd; ii++) {
	bet[ii].r -= beta[ii];
    }
    trl_print_complex_(info, title, jnd, bet, 1);
    strcpy(title, "Error in Lanczos recurrence (overall) =");
    trl_print_complex_(info, title, 1, aq, 1);
    //if( info->verbose > 7) {
    strcpy(title,
	   "|| A q_i - alpha_i q_i - beta_{i-1} q_{i-1} - beta_i q_{i+1} ||..");
    trl_print_complex_(info, title, jnd, wrk, 1);
    strcpy(title,
	   "(A q_i - alpha_i q_i - beta_{i-1} q_{i-1})*q_{i+1}/beta_i ..");
    trl_print_complex_(info, title, jnd, cs, 1);
    strcpy(title, "Sine of the angles ..");
    for (ii = 0; ii < jnd; ii++) {
	cs[ii].r = (cs[ii].r * cs[ii].r + cs[ii].i * cs[ii].i);
	if (cs[ii].r < one) {
	    cs[ii].r = sqrt(one - cs[ii].r);
	} else {
	    cs[ii].r = -one;
	}
    }
    trl_print_complex_(info, title, jnd, cs, 1);
    //}
    if (lwrk < jnd * 4 + nrow)
	Free(aq);
//
// .. end of trl_check_recurrence ..
//
}
