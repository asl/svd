/*
// ZTRLan routine 
// Lawrence Berkeley National Lab.
//
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
//
// ..
// .. internal libraries ..
#include "dsort2_i.h"
#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trlcore_i.h"
#include "trl_comm_i.h"

#include "ztrlcore_i.h"
#include "ztrlaux_i.h"
#include "ztrl_comm_i.h"
#include "ztrlan_i.h"

#define TUNED 0
// 0: no tuning is applied
// 1: checking for info->verbose is omitted
// 2: performance measurments are not taken
//
////
//
// Following are the internal subroutines for printing etc used in function lanczos_
//
////
void zlog_error_state(trl_info * info, int kept, int j1, int j2,
		      int jnd, int nrow, int mev, double *eval,
		      double *alpha, double *alfrot, double *beta,
		      double *betrot, trl_dcomplex * evec,
		      trl_dcomplex * base, trl_dcomplex * qa,
		      trl_dcomplex * qb, trl_dcomplex * rr, char *title,
		      int *iwrk)
{
// Purpose
// =======
// Dump important variables when return due to error
//
// Arguments
// =========
// info       (input) Pointer to structure trl_info_
//             On entry, points to the current TRL_INFO.
//
// kept       (input) Integer
//             On entry, specifies the number of lanczos vector kept at the last restart.
//
// j1         (input) Integer
//             On entry, specifies the last column of evec, that contains a lanczos vector.
//
// j2         (input) Integer
//             On entry, specifies the last column of base, that contains a lanczos vector.
//
// jnd        (input) Integer
//             On entry, specifies the number of lanczos vectors computed.
//
// nrow       (input) Integer
//             On entry, specifies the number of rows in the eigenvectors.
//
// mev        (input) Integer
//             On entry, specifies the maximum number of eigenvalues allowed.
//
// eval       (input) Double precision array of dimension (mev)
//             On entry, contains the eigenvalues computed.
//
// alpha      (input) Double precision array of dimension (info->maxlan)
//             On entry, contains the values of alpha computed.
//
// alfrot     (input) Double precision array of dimension (info->maxlan)
//             On entry, contains the values of alpha after rotation (diagonalization).
//
// beta       (input) Double precision array of dimension (info->maxlan)
//             On entry, contains the values of beta computed.
//
// betrot     (input) Double precisino array of dimension (info->maxlan)
//             On entry, contains the values of beta after rotation (diagonalization).
//
// evec       (input) Double precision array of dimension (nrow,mev)
//             On entry, contains the eigevectors computed.
//
// base       (input) Double precision array of dimension (nrow,nbas)
//             On entry, contains the lanczos vectors, that did not fit in evec.
//
// qa         (input) Double precision array of dimension (nrow)
//             On entry, contains the lanczos vector from the last iteration.
//
// qb         (input) Double precision array of dimension (nrow)
//             On entry, contains the lanczos vector from the two iterations ago.
//
// rr         (input) Double precision array of dimension (nrow)
//             On entry, contains the current lanczos vector being computed.
//
// title      (workspace) String length of (STRING_LEN)
//             On entry, provides a space to store the title of the information printed out.
//
// iwrk       (workspace) Integer array of dimension ()
//
// ..
// .. local variables ..
    FILE *fp = info->log_fp;
//
// ..
// .. executable statements ..
    trl_time_stamp(fp);
    strcpy(title, "Dumping the content of the variables on error..");
    iwrk[0] = info->stat;
    trl_print_int(info, title, 1, iwrk, 1);
    trl_terse_info(info, fp);
    fprintf(fp, "This Lanczos iteration started with %d vectors.\n", kept);
    fprintf(fp, "There are %d (%d, %d) Lanczos vectors currently.\n", jnd,
	    j1, j2);
    if (jnd != j1 + j2)
	jnd = j1 + j2;
    if (jnd < 0 || jnd > info->klan)
	jnd = 0;
    strcpy(title, "Content of eval ..");
    trl_print_real(info, title, mev, eval, 1);
    if (jnd > 0) {
	sprintf(title, "Alpha(1:%d) .. ", jnd);
	trl_print_real(info, title, jnd, alpha, 1);
	sprintf(title, " Beta(1:%d) .. ", jnd);
	trl_print_real(info, title, jnd, beta, 1);
	sprintf(title, "Alfrot(1:%d) .. ", jnd);
	trl_print_real(info, title, jnd, alfrot, 1);
	sprintf(title, "betrot(1:%d) .. ", jnd);
	trl_print_real(info, title, jnd, betrot, 1);
    }
    if (j1 > 0) {
	strcpy(title, "the First row of evec ..");
	trl_print_complex_(info, title, j1, evec, nrow);
	sprintf(title, "row %d of evec ..", nrow);
	trl_print_complex_(info, title, j1, &evec[nrow - 1], nrow);
    }
    if (j2 > 0) {
	strcpy(title, "the First row of base ..");
	trl_print_complex_(info, title, j2, base, nrow);
	sprintf(title, "row %d of base ..", nrow);
	trl_print_complex_(info, title, j2, &base[nrow - 1], nrow);
    }
    if (qb != NULL) {
	sprintf(title, "Content of qb (q_%d) ..", jnd - 1);
	trl_print_complex_(info, title, nrow, qb, 1);
    }
    if (qa != NULL) {
	sprintf(title, "Content of qa (q_%d) ..", jnd);
	trl_print_complex_(info, title, nrow, qa, 1);
    }
    if (rr != NULL) {
	sprintf(title, "Content of rr (residual == q_%d) ..", jnd + 1);
	trl_print_complex_(info, title, nrow, rr, 1);
    }
    if (info->my_pe == 0) {
      Rprintf("TRLanczos returned with error\n");
      Rprintf("Contents of most variables are dumped to log file %s.\n",
              info->log_file);
    }
//
//  .. end of print_error_state_ ..
//
}

////
void zprint_restart_state(trl_info * info, char *title, int nrow,
			  int mev, double *alpha, double *beta,
			  double *betrot, trl_dcomplex * evec, int lde,
			  double *yy, int kept, int locked, int *iwrk,
			  double *wrk2, int i2, int jml)
{
// Purpose
// =======
// Print the current solution status to the log file.
//
// Arguments
// =========
// info     (input) Pointer to structure trl_info_
//           On entry, points to the current TRL_INFO.
//
// title    (workspace) String of length (STRING_LEN)
//           On entry, provides space to store title to print out.
//
// nrow     (input) Integer
//           On entry, specifies the number of lanczos vector kept at the restart.
//
// mev      (input) Integer
//           On entry, specifies the maximum number of eigenvalues allowed.
//
// alpha    (input) Double precision array of dimension (info->maxlan)
//           On entry, contains the value of alphas computed so far.
//
// beta     (input) Double precision array of dimension (info->maxlan)
//           On entry, contains the value of beta computed so far.
//
// betrot   (input) Double precision array of dimension (info->maxlan)
//           On entry, contains the value of beta rotated.
//
// evec     (input) Double precision array of dimension (nrow,mev)
//           On entry, contains the eigenvectors computed.
//
// yy       (input) Double precision array of dimension (nrow,jml)
//           On entry, contains the litz vectors of the tridiagonal matrix
//           computed after the previous restart.
//
// kept     (input) Integer
//           On entry, specifies the number of lanczos vector kept at the restart.
//
// locked   (input) Integer
//           On entry, specifies the number of eigenvalues converged so far.
//
// iwrk     (workspace) Integer array of dimension (4*maxlan)
//           Integer workspace used to..
//
// wrk2     (workspace) Double precision array of dimension
//           Double precision workspace used to
//
// i2       (input) Integer
//           On entry, specifies
//
// jml      (input) Integer
//           On entry, specifies the number of litz vectors computed for
//           the current restart.
//
// ..
// .. local parameters ..
    long c__1 = 1;
//
// ..
// .. local scalars ..
    int i, j1, j2;
// ..
// .. executable statements ..
    iwrk[0] = kept + locked;
    iwrk[1] = locked + i2;
    strcpy(title, "Number of saved and locked Ritz pairs ..");
    trl_print_int(info, title, 2, iwrk, 1);
    if (info->verbose > 2) {
	if (iwrk[1] == 0) {
	    strcpy(title, "Ritz values saved (ascending ordered) ..");
	} else {
	    strcpy(title, "Ritz values saved (may not be ordered) ..");
	}
	trl_print_real(info, title, kept + locked, alpha, 1);
	strcpy(title, "Residual norms of the saved Ritz pairs ..");
	for (i = 0; i < (kept + locked); i++) {
	    betrot[i] = fabs(beta[i]);
	}
	trl_print_real(info, title, kept + locked, betrot, 1);
    }
    if (info->verbose > 7) {
	for (j1 = 0; j1 < imin2(kept, info->verbose); j1++) {
	    for (j2 = 0; j2 <= j1; j2++) {
		wrk2[j2] =
		    trl_ddot(jml, &yy[j2 * jml], c__1, &yy[j1 * jml],
			     c__1);
	    }
	    wrk2[j1] = wrk2[j1] - 1;
	    sprintf(title, "Orthogonality level of y(%d) ..", j1 + 1);
	    trl_print_real(info, title, j1 + 1, wrk2, 1);
	}
    }
    if (info->verbose > 10) {
	for (j1 = 0; imin2(kept, info->verbose); j1++) {
	    sprintf(title, "eigenvector %d of Q'AQ ..", j1);
	    trl_print_real(info, title, jml, &yy[(j1 - 1) * jml], 1);
	}
    }
    if (info->verbose > 10) {
	int j1n = imin2(nrow, info->verbose);
	for (j1 = 0; j1 < imin2(kept + locked, mev); j1++) {
	    sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
	    trl_print_complex_(info, title, j1n, &evec[j1 * lde], 1);
	}
    }
//
//  .. end of print_restart_state ..
//
}

////
void zprint_final_state(trl_info * info, char *title, int nrow, int mev,
			double *eval, int lde, double *beta,
			trl_dcomplex * evec, double *yy, int kept,
			int jml)
{
// Purpose
// =======
// print the final state
//
// Arguments
// =========
// info    (input) Pointer to structure trl_info_
//          On entry, points to the current TRL_INFO.
//
// title   (workspace) String of length (STRING_LEN)
//          On entry, provides space to store the title of the information to print out.
//
// nrow    (input) Integer
//          On entry, specifies the number of rows in the eigenvectors.
//
// mev     (input) Integer
//          On entry, specifies the maximum number of eigenvalues allowed.
//
// eval    (input) Double precision array of dimension (mev)
//          On entry, contains the eigenvalues computed.
//
// beta    (input) Double precision array of dimension (info->maxlan)
//          On entry, contains the value of beta computed.
//
// evec    (input) Double precision array of dimension (nrow,mev)
//          On entry, contains the eigenvectors computed.
//
// yy      (input) Double precision array of dimension (nrow,jml)
//          On entry, contains the litz vectors computed at the last restart.
//
// kept    (input) Integer
//          On entry, specifies the number of lanczos vectors kept at the last restart.
//
// jml     (input) Integer
//          On entry, specifies the number of new lanczos vectors computed at the last 
//          restart.
//
// ..
// .. local scalars ..
    int j1;
// ..
// .. executable statements ..
    strcpy(title, "Final eigenvalues  (in ascending order)..");
    trl_print_real(info, title, kept, eval, 1);
    if (info->verbose > 4) {
	strcpy(title, "Final residual norms..");
	trl_print_real(info, title, kept, beta, 1);
    }
    if (info->verbose > 8) {
	for (j1 = 0; j1 < imin2(kept, info->verbose); j1++) {
	    sprintf(title, "Eigenvector %d of Q''AQ ..", j1);
	    trl_print_real(info, title, jml, &yy[j1 * jml], 1);
	}
    }
    if (info->verbose > 10) {
	int j1n = imin2(nrow, info->verbose);
	for (j1 = 0; j1 < imin2(kept, mev); j1++) {
	    sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
	    trl_print_complex_(info, title, j1n, &evec[j1 * lde], 1);
	}
    }
//
// .. end of print_final_state ..
//
}

////
void zwrite_checkpoint(trl_info * info, char *title, int nrow,
			double *alpha, double *beta, trl_dcomplex * evec,
			trl_dcomplex * base, int lde, int j1n, int jnd,
			int ldb, int j2n)
{
//
// Purpose
// =======
// Output a check point.
//
// Arguments
// =========
// info     (input) Pointer to structure trl_info_
//           On entry, points to the current TRL_INFO.
//
// title    (input) String of length TITLE_LEN.
//           On entry, provides space to store the title of the information to print out.
//
// nrow     (input) Integer
//           On entry, specifies the number of rows in the eigenvectors.
//
// alpha    (input) Double precision array of length (info->maxlan)
//           On entry, contains the values of alpha computed.
//
// beta     (input) Double precision array of length (info->maxlan)
//           On entry, contains the value of beta computed.
//
// evec     (input) Double precision array of length (nrow,mev)
//           On entry, contains the eigenvectors computed.
//
// base     (input) Double precision array of length (ldb,nbase)
//           On entry, contains the lanczos vectors, that did not fit in evec.
//
// lde      (input) Integer
//           On entry, specifies the leading dimension of evec.
//
// j1n      (input) Integer
//           On entry, specifies the column index of evec, that stores the current lanczos 
//           vector.
//
// jnd      (input) Integer
//           On entry, specifies the number of lanczos vector computed so far.
//
// ldb      (input) Integer
//           On entry, specifies the leading dimension of base.
//
// j2n      (input) Integer
//           On entry, specifies the column index of base, the stores the current lanczos 
//           vector.
//
// ..
// .. local variables ..
    int ii;
    clock_t c1, c2;
//
// ..
// .. executable statements ..
    trl_pe_filename(138, title, info->cpfile, info->my_pe, info->npes);
    c1 = clock();
    ii = ztrl_write_checkpoint(title, nrow, alpha, beta, evec, lde, j1n,
				base, ldb, j2n);
    c2 = clock();
    if (c2 > c1) {
	info->clk_out = info->clk_out + (c2 - c1);
    } else {
	info->clk_out = info->clk_out + ((info->clk_max - c1) + c2);
    }
    info->wrds_out = info->wrds_out + jnd * (nrow + nrow + 2) + nrow + 2;
    info->stat = trl_sync_flag(info->mpicom, ii);
//
//  .. end of zwrite_checkpoint ..
//
}

////
void
ztrlanczos(ztrl_matprod op, trl_info * info, int nrow, int mev, double *eval,
           trl_dcomplex * evec, int lde, trl_dcomplex * wrk, int lwrk,
           double *dwrk, int ldwrk, void *lparam)
{
// Purpose
// =======
// The actual work routine of restarted Lanczos program for real
// symmetric eigenvalue problems
//
// user may directly invoke this sunroutine but she/he is responsible
// for input correct arguments as needed
//
// 1) info needs to be initialized
// 2) if info->nec>0, the first nec elements of eval and first nec
//    columns of evec must contain valid eigenpairs
// 3) workspace of currect size
//    eval(mev)
//    evec(lde, mev) (lde >= nrow, mev >= ned)
//    base(nrow, info->maxlan-mev+1) (nrow>=nrow, not used if mev>maxlan)
//    wrk(lwrk) minimum amount of memory required by TRLANCZOS is
//    maxlan*(maxlan+10)
// 4) if log files are to be written, the user need to open files on IO
//    unit log_io so that the log gile may be written correctly.
// 5) the user must set the timing variable info->clk_tot and
//    info->clk_max using system_clock function call in order for this
//    subroutine to track timing results correctly
//
// Algorithm
// =========
//  0. initialize input vector
//  1. do while (more eigenvalues to compute .and. more MATVEC allowed)
//  2.    first step
//     o   alpha(k+1) = dot_product(q_{k+1}, Aq_{k+1})
//     o   rr = A*q_{k+1}-alpha(k+1)*q_{k+1}-\sum_{i=1}^{k} beta(i)q_i
//     o   (re-orthogonalization)
//  3.    do j = k+2, m
//     o     rr = Aq_j
//     o     alpha(j) = dot_product(q_j, rr)
//     o     rr = rr - alpha(j)*q_j - beta(j-1)*q_{j-1}
//     o     (re-orthogonalization)
//        end do j = k+2, m
//  4.    restarting
//     o   call dstqrb to decompose the tridiagonal matrix
//     o   perform convergence test
//     o   determine what and how many Ritz pairs to save
//     o   compute the Ritz pairs to be saved
//     end do while
//
// The re-orthogonalization procedure is implemented in trl_orth.  it
// produces a normalized vector rr that is guaranteed to be orthogonal
// to the current basis.  An error will be reported if it can not
// achieve its goal.
//
// Arguments
// =========
// ops     (input) Functin pointer.
//          On entry, points to the function that performs the matrix-vector operation.
//          The operator that defines the eigenvalue problem is expected to have
//          the following interface
//          void op(nrow, ncol, xin, ldx, yout, ldy)
//             nrow  (input) Integer
//                    On entry, specifies the number of rows in xin and xout.
//             ncol  (input) Integer
//                    On entry, specifies the number of columns in xin and xout.
//             xin   (input) double precision array of dimension (ldx,ncol)
//                    On entry, specifies the vector/vectors for which the matrix-vector
//                    is performed.
//             ldx   (input) Integer
//                    On entry, specifies the leading dimension of xin
//             yout  (output) Double precision array of diimension (ldy,ncol)
//                    On exit, specifies the resulting vector/vectors.
//             ldy   (input) Integer
//                    On entry, specifies the leading dimension of yout.
//
// info    (input) Pointer to structure trl_info_
//          On entry, points to the current TRL_INFO.
//
// nrow    (input) Integer
//          On entry, specifies the number of rows in eigenvectors.
//
// mev     (input) Integer
//          On entry, specifies the number of columns allocated to store eigenvectors.
//
// eval    (output) Double array of dimension (mev)
//          On successful exit, contains the eigenvalues.
//
// evec    (output) Double array of dimension (lde,mev)
//          On successful exit, contains the eigenvectors.
//
// lde     (input) Integer
//          On entry, specifies the leading dimension of evec.
//
// base    (workspace) Double precision array of dimension (nrow,nbas)
//          Used to hold the lanczos vectors not fit in evec, i.e., nbas=info->maxlan-mev+1.
//
// nbas    (input) Integer
//          On entry, specifies the number of columns in base.
//
// wrk     (workspace) Double precision array of dimension (lwrk)
//          Workspace for lanczos iterations.
//
// lwrk    (input) Integer
//          On entry, specifies the size of workspace provided.
//
//
// ..
// .. local parameters ..
    int c__1 = 1;
//
// ..
// .. local variables ..
    char title[TRLAN_STRING_LEN];
    int i, i1, i2, j1, j2, jnd, jml, j1n, j2n, kept, prek, nbas, ldqa,
	ldrr, count;
    int next_test, chkpnt, locked;
    clock_t clk1;
    int *iwrk;
    int i__1;
    trl_dcomplex z__1;
    trl_dcomplex *rr, *qa, *qb, *base, *yy2;
    double *alpha, *beta, *alfrot, *betrot, *lambda, *res, *yy, *rot;
//
// ..
// .. executable statements ..
//
// initialize title and clock.
    strcpy(title, "");
    clk1 = 0;
//
// alpha, beta, alfrot and betrot have fixed locations in wrk, i.e.,
//   alpha: wrk(1:maxlan), beta: wrk(maxlan+1:2*maxlan),
//   alfrot: wrk(2*maxlan+1:3*maxlan), betrot: wrk(3*maxlan+1:4*maxlan)
//
    nbas = nrow * imax2(1, info->maxlan - mev + 1);
    base = Calloc(nbas, trl_dcomplex);
    alpha = Calloc(info->maxlan, double);
    beta = Calloc(info->maxlan, double);
    alfrot = Calloc(info->maxlan, double);
    betrot = Calloc(info->maxlan, double);
    lambda = Calloc(info->maxlan, double);
    res = Calloc(info->maxlan, double);
    yy = Calloc((info->maxlan) * (info->maxlan),
                double);
    yy2 = Calloc((info->maxlan) * (info->maxlan),
                 trl_dcomplex);
    rot = Calloc((info->maxlan) * (info->maxlan),
                 double);
    i2 = 0;
//
// allocate an integer workspace. iwrk holds...
    iwrk = Calloc((4 * info->maxlan), int);
//
// chkpnt specifies how often the check point should be written.
    if (info->cpflag <= 0) {
	// check point is not written
	chkpnt = info->maxmv + info->maxlan;
    } else {
	// check point is written at every chpnt matrix-vector operations.
	chkpnt = info->maxmv / info->cpflag;
    }
//
// locked specifies the number of eigenvalues converged.
    locked = info->nec;
//
// assign values to alpha, beta
// uses the assumption that the content of eval(1:nec) are eigenvalues
// and their residual norms are zero
//
    if (locked > 0) {
	memcpy(alpha, eval, locked * sizeof(double));
	memset(beta, 0, locked * sizeof(double));
    }
//
// get valid initial guess for the Lanczos iterations
// wrk2 points to the end of available wrk (first 4*maxlan hold alpha, beta, alfrot, and betrot)
// to the end of wrk.
// *** j1 and j2 are the number of colums, and not indices ***
    ztrl_initial_guess(nrow, evec, lde, mev, base, nrow, nbas, alpha,
                       beta, &j1, &j2, info, wrk, lwrk, lparam);
//
// On return from trl_initial_guess, j1 is the last column index of evec used, and
// j2 is the last column index of base used, i.e., jnd specifies the size of the current lanczos basis.
    jnd = j1 + j2;
    kept = jnd;
    if (info->stat != 0) {
	qa = NULL;
	qb = NULL;
	rr = NULL;
	kept = 0;
	j1 = 0;
	j2 = 0;
	jnd = 0;
	goto end;
    }
// we will perform the first convergence test after next_test
// matrix-vector multiplications
    i1 = info->ned - jnd;
    next_test = i1 + imin2(i1, imin2(6, info->ned / 2));
//
//*************************************************//
//            -- the TRLan outer loop --           //
// restart if                                      //
//    1. there is more eigenvalues to compute, and //
//    2. more matrix-vector operations are allowed //
//*************************************************//
//
    count = 0;
    while (info->matvec < info->maxmv && info->nec < info->ned) {
	//
	// jnd is the size of the current lanczos basis, so increment it.
	jnd++;
	//
	// the first iteration of TRLan
	//
	// qa points to the last lanczos base vector.
	if (j1 < mev) {
	    // there is still enough space in evec, so use it.
	    j1 = j1 + 1;
	    qa = &evec[(j1 - 1) * lde];
	    ldqa = lde;
	} else {
	    // no more space in evec, so use base.
	    j2 = j2 + 1;
	    qa = &base[(j2 - 1) * nrow];
	    ldqa = nrow;
	}
	//
	// j1n and j2n specify the location of the next lanczos basis, and
	// rr points the next lanczos base vector.
	if (j1 < mev) {
	    // there is still enough space in evec, so use it.
	    j1n = j1 + 1;
	    j2n = 0;
	    rr = &evec[(j1n - 1) * lde];
	    ldrr = lde;
	} else {
	    // no more space in evec, so use base.
	    j1n = mev;
	    j2n = j2 + 1;
	    rr = &base[(j2n - 1) * nrow];
	    ldrr = nrow;
	}
	//
	// perform matrix-vector multiplication, i.e., rr = A*qa
	// record the total time and the time in MATVEC

	clk1 = clock();
	if (clk1 < info->clk_tot) {
            info->tick_t +=
                (info->clk_max -
                 info->clk_tot) / (double) (info->clk_rate);
            info->tick_t +=
                (info->clk_max + clk1) / (double) (info->clk_rate);
            info->clk_tot = clk1;
	} else if (info->clk_tot < 0 && clk1 >= 0) {
	    info->tick_t -= info->clk_tot / (double) (info->clk_rate);
	    info->tick_t += clk1 / (double) (info->clk_rate);
	    info->clk_tot = clk1;
	}
	//
	i__1 = 1;
	op(&nrow, &i__1, qa, &ldqa, rr, &ldrr, lparam);

#if TUNED < 2
	add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
#endif

	(info->matvec)++;
	//
	// computed the next alpha = qa' * A * qa
	//printf( "computing alpha (%d)\n",nrow );
	trl_zdotc(&z__1, nrow, qa, c__1, rr, c__1);
	alpha[jnd - 1] = z__1.r;
	trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], dwrk);
	//
	// Perform the Lanczos orthogonalization.
	// rr = rr - sum_{i=1,...j1} beta(i)*evec(:,i) - sum_{1,...,j2} beta(j1+i)*base(:,i)
	// Just for a convenience beta[jnd-1]=alpha[jnd-1] just computed.
	//
	// orthogonalize with lanczos vectors stored in evec, first.
	beta[jnd - 1] = alpha[jnd - 1];
#if TUNED < 2
	info->flop = info->flop + nrow + nrow;
#endif

	if (j1 > 1) {
	    //
	    // compute rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(1),...,beta(i1)]'
	    zdgemv_(nrow, j1, -1.0, evec, lde, beta, 1.0, rr);

#if TUNED < 2
	    info->flop = info->flop + 2 * j1 * nrow;
#endif

	} else if (j1 == 1) {
	    // there is no beta, so just compute
	    //   rr = rr - alpha(1)*qa
	    zdaxpy_(nrow, -alpha[0], qa, rr);

#if TUNED < 2
	    info->flop = info->flop + nrow + nrow;
#endif
	}
	//
	// orthogonalize with lanczos vectors stored in base, now.
	if (j2 > 1) {
	    //
	    // compute rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(j1+1),...,beta(j1+j2)]'
	    zdgemv_(nrow, j2, -1.0, base, nrow, &beta[j1], 1.0, rr);

#if TUNED < 2
	    info->flop = info->flop + 2 * j2 * nrow;
#endif
	} else if (j2 == 1) {
	    // there is no beta, so just compute
	    zdaxpy_(nrow, -beta[jnd - 1], qa, rr);

#if TUNED < 2
	    info->flop = info->flop + nrow + nrow;
#endif
	}
	//
	// perform re-orthogonalization (full-orthogonalization)
#if TUNED < 2
	info->flop_h = info->flop_h - info->flop;
	clk1 = clock();
#endif

	ztrl_orth(nrow, evec, lde, j1, base, nrow, j2, rr, kept, alpha,
		   beta, wrk, lwrk, info);

#if TUNED < 1
	if (info->verbose > 8) {
	    // check orthogonality after the initilization step
	    ztrl_check_orth(info, nrow, evec, lde, j1n, base, nrow, j2n,
                      wrk, lwrk, lparam);
	}
#endif

#if TUNED < 2
	add_clock_ticks(info, &(info->clk_orth), &(info->tick_h), clk1);
	info->flop_h = info->flop_h + info->flop;
#endif

	if (info->stat != 0)
	    goto end;

#if TUNED < 1
	if (info->verbose > 5) {
	    print_alpha_beta(info, title, jnd, alpha, beta);
	}
#endif

	//
	// transform the matrix formed by alpha and beta into a
	// tridiagonal matrix, rot stores the transformation matrix
	//
	// the already-converged part is just diagonal.
	memcpy(alfrot, alpha, locked * sizeof(double));
	memset(betrot, 0, locked * sizeof(double));
	//
	// now, diagonalize the rest of matrix.
	i1 = jnd - locked;
	trl_tridiag(i1, &alpha[locked], &beta[locked], rot,
		     &alfrot[locked], &betrot[locked], dwrk, ldwrk,
		     &(info->stat));

#if TUNED < 2
	info->flop = info->flop + 8 * i1 * i1 * i1 / 3;	// Golub:1996:MC, P415
#endif

	if (info->stat != 0)
	    goto end;
	betrot[jnd - 1] = beta[jnd - 1];
	alfrot[jnd - 1] = alpha[jnd - 1];
	//
	//******************************************************//
	// regular iterations of Lanczos algorithm (inner loop) //
	// loop if                                              //
	//   1. there is space to store lanczos basis           //
	//   2. there is more eigenvalues to compute, and       //
	//******************************************************//
	//
	while (jnd < info->klan && info->nec < info->ned) {
	    //
	    // compute the kth lanczos vector.
	    //  qb is (k-2)nd lanczos vector, and qa is (k-1)st lanczos vector.
	    qb = qa;
	    qa = rr;
	    //
	    // increment j1, j2, and jnd.
	    j1 = j1n;
	    j2 = j2n;
	    jnd++;
	    //
	    // find the next available space for the kth lanczos vector.
	    if (j1n < mev) {
		// there is still a space in evec.
		j1n++;
		rr = &evec[(j1n - 1) * lde];
	    } else {
		// no more space in evec, so use base.
		j2n++;
		if (j2n <= nbas) {
		    rr = &base[(j2n - 1) * nrow];
		} else {
		    info->stat = -1111;
		    goto end;
		}
	    }
	    //
	    // perform the matrix-vector operation.

#if TUNED < 2
	    clk1 = clock();
#endif

	    i__1 = 1;
	    op(&nrow, &i__1, qa, &ldqa, rr, &ldrr, lparam);

#if TUNED < 2
	    add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
#endif

	    info->matvec = info->matvec + 1;
	    //
	    // compute alpha(jnd) = qa' * A * qa
	    trl_zdotc(&z__1, nrow, qa, c__1, rr, c__1);
	    alpha[jnd - 1] = z__1.r;
	    trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], dwrk);
	    //
	    // the Lanczos orthogonalization (three-term recurrence).
	    //   rr = rr - alpha(jnd)*qa - beta(jnd-1)*qb
	    zdaxpy_(nrow, -alpha[jnd - 1], qa, rr);
	    zdaxpy_(nrow, -beta[jnd - 2], qb, rr);

#if TUNED < 2
	    info->flop = info->flop + 6 * nrow;
#endif

	    //
	    // re-orthogonalization, and compute beta(jnd)

#if TUNED < 2
	    info->flop_h = info->flop_h - info->flop;
	    clk1 = clock();
#endif

	    ztrl_orth(nrow, evec, lde, j1, base, nrow, j2, rr, kept,
		       alpha, beta, wrk, lwrk, info);

#if TUNED < 2
	    add_clock_ticks(info, &(info->clk_orth), &(info->tick_h),
			     clk1);
	    info->flop_h = info->flop_h + info->flop;
#endif
	    //
	    // copy alpha and beta into alfrot and betrot
	    alfrot[jnd - 1] = alpha[jnd - 1];
	    betrot[jnd - 1] = beta[jnd - 1];
	    if (info->stat != 0)
		goto end;

#if TUNED < 1
	    if (info->verbose > 4) {
		print_alpha_beta(info, title, jnd, alpha, beta);
	    }
#endif

	    //
	    // perform convergence test once in a while
	    if (info->matvec >= next_test) {

#if TUNED < 1
		if (info->verbose > 5) {
		    print_all_alpha_beta(info, title, jnd, alfrot,
					  betrot);
		}
#endif

		trl_get_eval(jnd, locked, alfrot, betrot, lambda, res,
			      dwrk, ldwrk, &(info->stat));
		if (info->stat != 0)
		    goto end;

#if TUNED < 1
		if (info->verbose > 2) {
		    print_lambda_res(info, jnd, lambda, res);
		}
#endif

		i1 = imin2(mev, jnd);
		memcpy(eval, dwrk, i1 * sizeof(double));
		trl_convergence_test(jnd, lambda, res, info, dwrk);
		//
		// decide when to perform the next test
		if (info->nec < info->ned && info->nec > 0) {
		    //
		    // assuming a same number of matrix-vector product is required for each eigenvalues to converge.
		    next_test =
			(double) (info->ned * info->matvec) /
			(double) (info->nec);
		} else if (info->nec == 0) {
		    next_test = next_test + next_test;
		    if (info->maxlan == info->ntot) {
			next_test =
			    (int) ceil(0.5 *
				       (info->maxlan + info->matvec));
		    }
		}
#if TUNED < 1
		if (info->verbose > 0)
		    trl_print_progress(info);
#endif
	    }
	}
	//clk2 = clock();
	//***************************************************************//
	// end of inner (regular Lanczos three-term recurrence) loop     //
	//***************************************************************//
	//
	// error checking for debugging use
	//

#if TUNED < 1
	if (info->verbose > 6) {
	    ztrl_check_orth(info, nrow, evec, lde, j1n, base, nrow, j2n,
                      wrk, lwrk, lparam);
	    if (info->verbose > 7) {
		ztrl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
				       base, nrow, j2n, kept, alpha, beta,
                          wrk, lwrk, lparam);
	    }
	}
#endif
	//
	// convert the integer counters to floating-point counters
	//

#if TUNED < 2
	i2 = info->clk_max / 4;
	if (info->flop > i2) {
	    info->rflp = info->rflp + info->flop;
	    info->flop = 0;
	}
	if (info->flop_h > i2) {
	    info->rflp_h = info->rflp_h + info->flop_h;
	    info->flop_h = 0;
	}
	if (info->flop_r > i2) {
	    info->rflp_r = info->rflp_r + info->flop_r;
	    info->flop_r = 0;
	}
	if (info->clk_op > i2) {
	    info->tick_o =
		info->tick_o + info->clk_op / (double) (info->clk_rate);
	    info->clk_op = 0;
	}
	if (info->clk_orth > i2) {
	    info->tick_h =
		info->tick_h + info->clk_orth / (double) (info->clk_rate);
	    info->clk_orth = 0;
	}
	if (info->clk_res > i2) {
	    info->tick_r =
		info->tick_r + info->clk_res / (double) (info->clk_rate);
	    info->clk_res = 0;
	}
	info->flop_r = info->flop_r - info->flop;
	//
	// *** Determine whether to restart ***
	// compute the Ritz values and Ritz vectors if they are not up to
	// date
	//
	clk1 = clock();
#endif

	prek = kept;
	jml = jnd - locked;
	i2 = kept - locked + 1;
	if (info->nec < info->ned) {
	    // need to compute the updated Ritz values and residual norms
	    //
	    // Given tridiagonal matrix (diagonals stored in alfrot, and off-diagonals stored in betrot),
	    // computes Ritz value (approximate eigenvalues), using dstqrb, and retrned in lambda.
	    // the last components of the eigenvectors are also computed.
	    //clk3 = clock();
	    trl_get_eval(jnd, locked, alfrot, betrot, lambda, res, dwrk,
			  ldwrk, &(info->stat));
	    if (info->stat != 0)
		goto end;

#if TUNED < 1
	    if (info->verbose > 2) {
		print_lambda_res(info, jnd, lambda, res);
	    }
#endif

	    //clk3 = clock();
	    trl_convergence_test(jnd, lambda, res, info, dwrk);
	    //printf( "conv test: %d\n",(int)(clock()-clk3) );

#if TUNED < 1
	    if (info->verbose > 0) {
		trl_print_progress(info);
	    }
#endif
	}
	//
	// Given the tridiagonal matrix and Ritz values, compute the Ritz vectors
	// (rotational matrix, used for tridiagonalization, is also applied).
	// Also, decide how many vectors to save if restart
	if (info->nec < info->ned && info->matvec < info->maxmv) {
	    //
	    // prepare to restart, reorder the eigenvales based on the input parameter.
	    //clk3 = clock();
	    trl_shuffle_eig(jml, info->klan - locked, &lambda[locked],
			     &res[locked], info, &kept, locked);
	    //
	    // compute eigenvectors using dstein (inverse interations)
	    //clk3 = clock();
	    if (kept * 3 < jml) {
		trl_get_tvec(jml, &alfrot[locked], &betrot[locked], 0, i2,
			      rot, kept, &lambda[locked], yy, iwrk, dwrk,
			      ldwrk, &(info->stat));
		if (info->stat == 0 && (locked + kept) > 0) {
		    memcpy(alpha, lambda,
			   (locked + kept) * sizeof(double));
		}
	    }
	    //
	    // compute eigenvectors using dsyev (QR)
	    //clk3 = clock();
	    if (kept * 3 >= jml || info->stat != 0) {
		if ((locked + kept) > 0)
		    memcpy(alfrot, lambda,
			   (locked + kept) * sizeof(double));
		trl_get_tvec_a(jml, prek - locked, &alpha[locked],
				&beta[locked], kept, &alfrot[locked], yy,
				dwrk, ldwrk, iwrk, &(info->stat));
	    }
	    if (info->stat != 0)
		goto end;
	    for (i = 0; i < kept; i++) {
		beta[locked + i] = yy[(i + 1) * jml - 1] * betrot[jnd - 1];
	    }
	    if (jml > info->ned + (info->ned / 5 + 6)) {
		trl_set_locking(jml, kept, &alpha[locked], &beta[locked],
				 yy, info->anrm, &i2);
	    } else {
		i2 = 0;
	    }
	    //
	    // generate Ritz vectos, reclaim the space pointed by ROT
	    //clk3 = clock();
	    ztrl_ritz_vectors(nrow, locked, kept, yy, jml, evec, lde, j1,
			       base, nrow, j2, wrk, lwrk, yy2);

#if TUNED < 2
	    info->flop = info->flop + 2 * nrow * jml * kept;
#endif

#if TUNED < 1
	    if (info->verbose > 0) {
		zprint_restart_state(info, title, nrow, mev, alpha, beta,
				      betrot, evec, lde, yy, kept, locked,
				      iwrk, dwrk, ldwrk, jml);
	    }
#endif
	    //
	    // reset the counters and indices to the correct values for restarting
	    kept += locked;
	    locked += i2;
	    info->locked = locked;
	    jnd = kept;
	    if (jnd <= mev) {
		j1 = jnd;
		j2 = 0;
	    } else {
		j1 = mev;
		j2 = jnd - mev;
		if (j2 >= (nbas - 1)) {
		    info->stat = -1111;
		    goto end;
		}
	    }
	    if (info->nec > 0) {
		next_test =
		    (int) (((double) (info->matvec * info->ned)) /
			   ((double) info->nec));
	    } else {
		next_test = next_test + info->maxlan;
	    }
	    i1 = imin2(mev, jnd);
	    if (i1 > 0)
		memcpy(eval, lambda, i1 * sizeof(double));
	    if (jnd < mev) {
		j1n = j1 + 1;
		j2n = 0;
		//memcpy( &evec[(j1n-1)*nrow], rr, nrow*sizeof(trl_dcomplex) );
		memcpy(&evec[(j1n - 1) * lde], rr,
		       nrow * sizeof(trl_dcomplex));
	    } else {
		j1n = mev;
		j2n = j2 + 1;
		memcpy(&base[(j2n - 1) * nrow], rr,
		       nrow * sizeof(trl_dcomplex));
	    }
	    //
	    // write checkpoint files
	    if (info->matvec >= chkpnt) {
		zwrite_checkpoint(info, title, nrow, alpha, beta, evec,
				   base, lde, j1n, jnd, nrow, j2n);
		chkpnt = chkpnt + info->maxmv / info->cpflag;
	    }
	} else {
	    //
	    // all wanted eigenpairs converged or maximum MATVEC used
	    // sort the eigenvalues in final output order
	    kept = imin2(info->nec, imax2(info->ned, mev - 1));
	    info->nec = kept;
	    if (kept == 0)
		kept = imin2(mev - 1, info->ned);
	    trl_sort_eig(jnd, info->lohi, kept, info->ref, lambda, res);
	    memcpy(eval, lambda, kept * sizeof(double));
	    if (kept * 3 < jnd) {
		// eigenvectors of the projection matrix (try inverse
		// interations)
		trl_get_tvec(jnd, alfrot, betrot, locked, i2, rot,
			      kept, eval, yy, iwrk, dwrk, ldwrk,
			      &(info->stat));
	    }
	    if (kept * 3 >= jnd || info->stat != 0) {
		//
		// too many eigenvectors or inverse iterations have failed,
		// try QR
		trl_get_tvec_a(jnd, prek, alpha, beta, kept, eval, yy,
				dwrk, ldwrk, iwrk, &(info->stat));
		if (info->stat != 0)
		    goto end;
	    }
	    if (kept > 0)
		memcpy(alpha, eval, kept * sizeof(double));
	    for (i = 0; i < kept; i++) {
		beta[i] = betrot[jnd - 1] * yy[(1 + i) * jnd - 1];
	    }
	    //
	    // generate eigenvectos, reclaim the space pointed by ROT
	    ztrl_ritz_vectors(nrow, 0, kept, yy, jnd, evec, lde, j1, base,
			       nrow, j2, wrk, lwrk, yy2);

#if TUNED < 2
	    info->flop = info->flop + 2 * nrow * jml * kept;
#endif

#if TUNED < 1
	    if (info->verbose > 1) {
		zprint_final_state(info, title, nrow, mev, eval, lde,
				    beta, evec, yy, kept, jml);
	    }
#endif

	    // reset the counters and indices to be used by check_orth and
	    // check_recurrence
	    jnd = kept;
	    j1 = kept;
	    j2 = 0;
	    if (j1 < mev) {
		j1n = j1 + 1;
		j2n = 0;
		//memcpy( &evec[(j1n-1)*nrow], rr, nrow*sizeof(trl_dcomplex) );
		memcpy(&evec[(j1n - 1) * lde], rr,
		       nrow * sizeof(trl_dcomplex));
	    } else {
		j1n = mev;
		j2n = 1;
		memcpy(base, rr, nrow * sizeof(trl_dcomplex));
	    }
	    // write checkpoint files
	    if (info->cpflag > 0) {
		zwrite_checkpoint(info, title, nrow, alpha, beta, evec,
				   base, lde, j1n, jnd, nrow, j2n);
	    }
	}
	//
	// check the orthogonality of the basis vectors before restarting
	//

#if TUNED < 1
	if (info->verbose > 6) {
	    ztrl_check_orth(info, nrow, evec, lde, j1n, base, nrow, j2n,
                      wrk, lwrk, lparam);
	    if (info->verbose > 7) {
		ztrl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
				       base, nrow, j2n, kept, alpha, beta,
                          wrk, lwrk, lparam);
	    }
	}
#endif
	//printf( "total: %d\n",(int)(clock()-clk2) );

#if TUNED < 2
	add_clock_ticks(info, &(info->clk_res), &(info->tick_r), clk1);
	info->flop_r = info->flop_r + info->flop;
#endif
	info->nloop = info->nloop + 1;
    }

    //*********************//
    // end of restart_loop //
    //*********************//
    // write the estimated residual norms to the beginning of WRK
    for (i = 0; i < j1; i++) {
	wrk[i].r = fabs(beta[i]);
    }
  end:
    if (info->stat < 0 && (info->verbose > 0 || info->my_pe == 0)) {
	zlog_error_state(info, kept, j1, j2, jnd, nrow, mev, eval, alpha,
			  alfrot, beta, betrot, evec, base, qa, qb, rr,
			  title, iwrk);
    }

    Free(iwrk);
    Free(rot);
    Free(yy2);
    Free(yy);
    Free(res);
    Free(lambda);
    Free(betrot);
    Free(alfrot);
    Free(beta);
    Free(alpha);
    Free(base);

    return;
//
// .. end of lanczos_ ..
//
}

////
void trl_smooth_zz(int n, trl_dcomplex * rr)
{
//
// Purpose
// =======
// Smooth out a vector, i.e.,
//  rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1) in Fortran.
// Used in trl_initial_guess.
//
// Arguments
// =========
// n   (input) Integer
//      On entry, specifies the dimension of rr.
//
// rr  (input/output) Double precision array of dimension (n)
//      On entry, contains the initial state of the vector. On exit, contain the vector after
//      the smothing is applied.
//
// ..
// .. local scalars ..
    int i;
// ..
// .. executable statements ..
    if (n <= 0)
	return;
    double r1, i1, r2, i2;
    r2 = rr[0].r;
    i2 = i1 = rr[0].i;
    rr[0].r = 2 * rr[0].r + rr[2].r + rr[n - 1].r;
    rr[0].i = 2 * rr[0].i + rr[2].i + rr[n - 1].i;
    for (i = 1; i < n - 1; i++) {
	r1 = rr[i].r;
	i1 = rr[i].i;
	rr[i].r = 2 * rr[i].r + rr[i + 1].r + r2;
	rr[i].i = 2 * rr[i].i + rr[i + 1].i + i2;
	r2 = r1;
	i2 = i1;
    }
    rr[n - 1].r = 2 * rr[n - 1].r + rr[1].r + r2;
    rr[n - 1].i = 2 * rr[n - 1].r + rr[1].i + i2;
//
// .. end of trl_smooth_zz ..
//
}

////
void ztrl_initial_guess(int nrow, trl_dcomplex * evec, int lde, int mev,
                        trl_dcomplex * base, int ldb, int nbas,
                        double *alpha, double *beta, int *j1, int *j2,
                        trl_info * info, trl_dcomplex * wrk, int lwrk,
                        void *lparam)
{
//
// Purpose
// =======
// check to make sure the initial guess vector contains valid nonzero numbers if not fill with 
// random numbers this routine will also read the checkpoint files to restore the previous state 
// of the Lancozs iterations
//
// Arguments
// =========
// nrow   (input) integer
//         On entry, specifies the number of rows in eigenvectors.
//
// evec   (input/output) double array of dimension (lde,mev)
//         On entry, the (nec+1)st column contains the initial guess.
//
// lde    (input) integer
//         On entry, specifies the leading dimention of evec.
//
// mev    (input) integer
//         On entry, specifies the number of Ritz vectors that can be stored in evec.
//
// base   (input/output) double array of dimension (ldb,nbas)
//         Stores the Ritz vectors, that cannot be stored in evec.
//
// ldb    (input) integer
//         On entry, specifies the leading dimention of base.
//
// nbas   (input) integer
//         On entry, specifies the number of Ritz vectors that can be stored in base
//
// alpha  (input/output) double array of dimension (mev+nbas-1)
//         On exit, stores alpha values if checkpoint files are provided.
//
// beta   (input/output) double array of dimension (mev+nbas-1)
//         On exit, stores beta values if checkpoint files are provided.
//
// j1     (output) pointer to integer
//         On exit, stores j1 (number of Ritz vectors in evec) if checkpoint files are 
//         provided.
//
// j2     (output) pointer to integer
//         On exit, stores j1 (number of Ritz vectors in base) if checkpoint files are 
//         provided.
//
// info   (input/output) pointer to trl_info structure
//         On entry, points to the data structure to store the current information about
//         the eigenvalue problem and the progress of TRLAN.
//
// wrk    (workspace) double array of dimension (lwrk)
//
// lwrk   (input) integer
//         Specifies the size of the workspace.
//
//
// parameter
    long c__1 = 1;
    //
    // local variable
    //
    int i, j, k, nran, north;
    double tmp, rnrm;
    clock_t ii, jj;
    char file[TRLAN_STRING_LEN];
    //
    // generate random seeds based on current clock ticks
    ii = clock();

    if (info->my_pe > 0) {
	ii = ii - (int) (info->my_pe * sqrt((double) ii));
    }
    //
    j = info->nec;
    if (info->guess > 1) {
	// retrieve a check-point file
	i = info->cpio;
	if (info->oldcpf != 0 && strlen(info->oldcpf) > 0) {
	    trl_pe_filename(TRLAN_STRING_LEN, file, info->oldcpf, info->my_pe,
			     info->npes);
	} else {
	    trl_pe_filename(TRLAN_STRING_LEN, file, info->cpfile, info->my_pe,
			     info->npes);
	}

#if TUNED < 2
	ii = clock();
#endif

	//i = ztrl_read_checkpoint( file, nrow, &evec[j*nrow], lde, mev-info->nec, j1,
	//          base, ldb, nbas, j2, (mev+nbas-1-j), &alpha[j], (mev+nbas-1-j),&beta[j] );
	i = ztrl_read_checkpoint(file, nrow, &evec[j * lde], lde,
				  mev - info->nec, j1, base, ldb, nbas, j2,
				  (mev + nbas - 1 - j), &alpha[j],
				  (mev + nbas - 1 - j), &beta[j]);
	info->stat = trl_sync_flag(info->mpicom, i);

#if TUNED < 2
	jj = clock();
	if (jj > ii) {
	    info->clk_in = jj - ii;
	} else {
	    info->clk_in = (info->clk_max - ii) + jj;
	}
#endif

	info->wrds_in = (*j1 + *j2) * (nrow + nrow + 2) + nrow + 2;
	*j1 = *j1 + info->nec;
	if (info->stat != 0)
	    return;
    } else {
	if (info->guess <= 0) {
	    // generate an arbitrary initial starting vector
	    // if (info->guess == 0), use the vector [1, 1, ...]^T
	    // else perturb some random locations of the above vector
	    if (info->guess == 0) {
		for (k = 0; k < nrow; k++) {
		    evec[j * lde + k].r = 1.0;
		    evec[j * lde + k].i = 0.0;
		}
	    } else {
		for (k = 0; k < nrow; k++) {
		    evec[j * lde + k].r = 1.0;
		    evec[j * lde + k].i = 1.0;
		}
		nran = imin2(-info->guess, lwrk);
		nran = 2 * (nran / 2);
    GetRNGstate();
		if (nran < nrow) {
		    for (k = 0; k < nran; k++) {
          wrk[k].r = unif_rand();
          wrk[k].i = unif_rand();
		    }
		    for (i = 0; i < nran - 1; i += 2) {
          ii = (int) (nrow * wrk[i].r);
          evec[j * lde + ii].r =
            evec[j * lde + ii].r + wrk[i + 1].r - 0.5;
          evec[j * lde + ii].i =
            evec[j * lde + ii].i + wrk[i + 1].i - 0.5;
		    }
		    info->flop += 4 * nran;
		} else {
      for (i = 0; i < nrow; i++) {
        evec[j * lde + i].r = unif_rand();
        evec[j * lde + i].i = unif_rand();
      }
      //trl_smooth_zz( nrow,&evec[(info->nec)*nrow] );
      trl_smooth_zz(nrow, &evec[(info->nec) * lde]);
      info->nrand++;
      info->flop += 8 * nrow;
		}
    PutRNGstate();
	    }
	}
	*j1 = info->nec;
	*j2 = 0;
    }
    tmp = 0.0;
    // make sure the norm of the next vector can be computed
    //zdotc_( wrk,&nrow,&evec[j*nrow],&c__1,&evec[j*nrow],&c__1 );
    trl_zdotc(wrk, nrow, &evec[j * lde], c__1, &evec[j * lde], c__1);
    trl_g_sum(info->mpicom, 1, &(wrk[0].r), &(wrk[1].r));

#if TUNED < 2
    info->flop = info->flop + nrow + nrow;
#endif

    if (wrk[0].r >= DBL_MIN && wrk[0].r <= DBL_MAX) {
	// set rnrm to let trl_CGS normalize evec(1:nrow, j)
	rnrm = sqrt(wrk[0].r);
    } else {
      GetRNGstate();
      for (i = 0; i < nrow; i++) {
        evec[j * lde + i].r = unif_rand();
        evec[j * lde + i].i = unif_rand();
      }
      PutRNGstate();
  
	//trl_smooth_zz( nrow,&evec[(info->nec)*nrow] );
	trl_smooth_zz(nrow, &evec[(info->nec) * lde]);
	info->nrand++;

#if TUNED < 2
	info->flop += 8 * nrow;
#endif
    }
    //
    // orthogonalize initial guess against all existing vectors
    //
    i = 0;
    tmp = 1.0;

#if TUNED < 2
    nran = info->nrand;
    north = info->north;
#endif

    if (*j1 < mev) {
	info->stat = ztrl_cgs(info, nrow, evec, lde, *j1, base, ldb, 0,
			       &evec[(*j1) * lde], &rnrm, &tmp, &i, wrk);

    } else if (*j2 <= 0) {
	info->stat = ztrl_cgs(info, nrow, evec, lde, *j1, evec, lde, 0,
			       base, &rnrm, &tmp, &i, wrk);

    } else {
	info->stat = ztrl_cgs(info, nrow, evec, lde, *j1, base, ldb, *j2,
			       &base[(*j2) * nrow], &rnrm, &tmp, &i, wrk);
    }

#if TUNED < 2
    info->flop =
	info->flop + 8 * nrow * ((info->north - north) * (*j1 + *j2) +
				 info->nrand - nran) + nrow;
#endif

#if TUNED < 1
    if (info->verbose > 6) {
	if (*j1 < mev) {
	    i = *j1 + 1;
	    ii = *j2;
	} else {
	    i = *j1;
	    ii = *j2 + 1;
	}
	ztrl_check_orth(info, nrow, evec, lde, *j1, base, ldb, ii, wrk,
                  lwrk, lparam);
    }
#endif
    return;
//
// .. end of trl_initial_guess ..
//
}

////
void ztrl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
			trl_dcomplex * vec1, int ld1, int m1,
			trl_dcomplex * vec2, int ld2, int m2,
			trl_dcomplex * wrk, int lwrk, trl_dcomplex * yy2)
{
//
// Purpose
// =======
// compute the Ritz vectors from the basis vectors and the eigenvectors of the projected system
// the basis vectors may be stored in two separete arrays the result need to be stored back in them
// lwrk should be no less than ny (lwrk>=ny) ***NOT checked inside***
//
// Arguments
// =========
// nrow   (input) Integer
//         On entry, specifies the number of rows in eigenvectors.
//
// lck    (input) Integer
//         On entry, specifies the number of Ritz values locked.
//
// ny     (input) Integer
//         On entry, specifies the number of columns in yy.
//
// yy     (input) Double precision array (ldy,ny)
//         On entry, contains the eigenvector of the "tri-diagonal" matrix.        
// 
// ldy    (input) Integer
//         On entry. specify the leading dimention of yy.
//
// vec1   (input) Complex array (ld1,m1)
//         On entry, contains the first part of Lanczos basis.
//
// m1     (input) Integer
//         On entry, specifies the number of Lanczos basis stored in vec1.
//
// ld1    (input) Integer
//         On entry, specifies the leading dimention of vec1.
//
// vec2   (input) Complex array (ld2,m2)
//         On entry, contains the second part of Lanczos basis.
//
// m2     (input) Integer
//         On entry, specifies the number of Lanczos basis stored in vec2.
//
// ld2    (input) Integer
//         On entry, specifies the leading dimention of vec2.
//
// wrk    (workspace) Complex array (lwrk)
// yy2    (workspace) Complex array (ldy,ny)
// lwrk   (input)
//         Specifies the size of the workspace.         
//
// ..
// .. local parameters ..
    int c__1 = 1;
    char notrans = 'N';
    trl_dcomplex zero, one;
//
// .. local variables ..
    int i, j, k, stride, ii, jl1, jl2, il1, il2, kv1, lwrk1;
//
// ..
// .. executable statements ..
    zero.r = 0.0;
    zero.i = 0.0;
    one.r = 1.0;
    one.i = 0.0;
    if (lck <= m1) {
	il1 = lck + 1;
	jl1 = m1 - lck;
	il2 = 1;
	jl2 = m2;
    } else {
	il1 = m1 + 1;
	jl1 = 0;
	il2 = lck - m1 + 1;
	jl2 = m1 + m2 - lck;
    }

    if (jl1 == 0 && jl2 == 0)
	return;

    lwrk1 = lwrk;
    for (i = 0; i < (ny * ldy); i++) {
	yy2[i].r = yy[i];
	yy2[i].i = 0.0;
    }

    kv1 = imin2(m1 - il1 + 1, ny);
    memset(wrk, 0, lwrk1 * sizeof(trl_dcomplex));
    //printf( "ny=%d\r\n",ny );
    if (ny > 1) {
	stride = lwrk1 / ny;
	for (i = 0; i < nrow; i += stride) {
	    j = imin2(nrow - 1, i + stride - 1);
	    k = j - i + 1;

	    //printf( "j1l=%d",jl1 );
	    if (jl1 > 0) {
		// compute wrk = vec1(i:j,:)*yy
		// (Note the leading dimension of vec1 is ld1. This effectively shift the vec1(i:j) to the top of the matrix.

		trl_zgemm(&notrans, &notrans, k, ny, jl1, one,
			  &vec1[(il1 - 1) * ld1 + i], ld1, yy2, ldy, zero,
			  wrk, k);

		/*
		   zdgemm_( k, ny, jl1, &vec1[(il1-1)*ld1+i], ld1, yy, ldy, wrk, k );
		 */
	    } else {
		memset(wrk, 0, lwrk1 * sizeof(trl_dcomplex));
	    }

	    //printf( "j2l=%d",jl2 );
	    if (jl2 > 0) {
		trl_zgemm(&notrans, &notrans, k, ny, jl2, one,
			  &vec2[(il2 - 1) * ld2 + i], ld2, &yy2[jl1], ldy,
			  one, wrk, k);

		/*
		   zdgemm2_( k, ny, jl2, &vec2[(il2-1)*ld2+i],ld2, &yy[jl1], ldy, wrk, k );
		 */
	    }

	    for (ii = 0; ii <= kv1 - 1; ii++) {
		memcpy(&vec1[(ii + il1 - 1) * ld1 + i], &wrk[ii * k],
		       k * sizeof(trl_dcomplex));
		//printf( "vec1=%e+%e*i\r\n",vec1[(ii+il1-1)*ld1+i].r,vec1[(ii+il1-1)*ld1+i].i );
	    }

	    for (ii = 0; ii <= (ny - kv1 - 1); ii++) {
		memcpy(&vec2[(ii + il2 - 1) * ld2 + i],
		       &wrk[(kv1 + ii) * k], k * sizeof(trl_dcomplex));
	    }
	}
    } else if (ny == 1) {
	stride = lwrk1;
	for (i = 0; i < nrow; i += stride) {
	    j = imin2(nrow - 1, i + stride - 1);
	    k = j - i + 1;
	    if (jl1 > 0) {
		/*
		   zdgemv_( k, jl1, 1.0, &vec1[(il1-1)*ld1+i], ld1, yy, 0.0, wrk );
		 */

		trl_zgemv(&notrans, k, jl1, one,
			  &vec1[(il1 - 1) * ld1 + i], ld1, yy2, c__1, one,
			  wrk, c__1);

		if (jl2 > 0) {

		    trl_zgemv(&notrans, k, jl2, one,
			      &vec2[(il2 - 1) * ld2 + i], ld2, &yy2[jl1],
			      c__1, one, wrk, c__1);

		    /*
		       zdgemv_( k, jl2, 1.0, &vec2[(il2-1)*ld2+i], ld2, &yy[jl1], 1.0 , wrk );
		     */
		}
	    } else {

		trl_zgemv(&notrans, k, jl2, one,
			  &vec2[(il2 - 1) * ld2 + i], ld2, &yy2[jl1], c__1,
			  zero, wrk, c__1);

		/*
		   zdgemv_( k, jl2, 1.0, &vec2[(il2-1)*ld2+i], ld2, &yy[jl1], 0.0, wrk );
		 */
	    }
	    if (kv1 > 0) {
		memcpy(&vec1[(il1 - 1) * ld1 + i], wrk,
		       k * sizeof(trl_dcomplex));
	    } else {
		memcpy(&vec2[(il2 - 1) * ld2 + i], wrk,
		       k * sizeof(trl_dcomplex));
	    }
	}
    }
//
// .. end of trl_ritzs_vectors_ ..
//
}

////
void ztrl_orth(int nrow, trl_dcomplex * v1, int ld1, int m1,
		trl_dcomplex * v2, int ld2, int m2, trl_dcomplex * rr,
		int kept, double *alpha, double *beta, trl_dcomplex * wrk,
		int lwrk, trl_info * info)
{
//
// Purpose
// =======
// Applies full re-orthogonalization;
//  1. if (global re-orthogonalization is needed)
//      call trl_cgs
//    else
//      perform extended local re-reorthogonalization
//    endif
//  2. perform normalization
//
// Arguments:
// ==========
// nrow   (input) Integer
//         On entry, specifies the number of rows in eigenvectors.
//
// v1     (input) Complex array (ld1,m1)
//         On entry, contains the first part of Lanczos basis computed.
//
// ld1    (input) Integer
//         On entry, specifies the leading dimention of v1.
//
// m1     (input) Integer
//         On entry, specifies the number of Lanczos basis in v1.
//
// v2     (input) Complex array (ld2,m2)
//         On entry, contains the second part of Lanczos basis computed.
//
// ld2    (input) Integer
//         On entry, specifies the leading dimention of v2.
//
// m2     (input) Integer
//         On entry, specifies the number of Lanczos basis in v2.
//
// rr     (input/output) Complex array (nrow)
//         On entry, contains the new Lanczos basis computed.
//         On exit, contains the next Lanczos basis computed after the orthogonalization.
//
// kept   (input) Integer
//         On etnry, specifies the number of Ritz vectors kept.
//
// alpha  (input/output) double precision array (m1+m2)
//         On entry, contains the alpha values, on exit, they are updated. 
//
// beta   (input/output) double precision array (m1+m2)
//         On entry, contains the beta values, on exit, they are updated if necessary,
//         (full orthogonalization).
//
// wrk    (workspace) complex array (lwrk)
//
// lwrk   (input) Integer
//         Specifies the size of workspace.
//
//
// ..
// .. local parameters ..
    double zero = 0.0, one = 1.0;
    long c__1 = 1;
//
// ..
// .. local variables ..
    trl_dcomplex z__1;
    int i, jnd, jm1, no, nr, usecgs;
    double tmp1, d__1;
    trl_dcomplex *qa, *qb;
//
// ..
// .. executable statements ..
//
// check for workspace size
    jnd = m1 + m2;
    jm1 = jnd - 1;
    if (ld1 >= nrow && ld2 >= nrow && lwrk >= imax2(4, jnd + jnd)) {
	info->stat = 0;
    } else {
	info->stat = -101;
	return;
    }
//
// compute the norm of the vector RR
//
    trl_zdotc(&(wrk[0]), nrow, rr, c__1, rr, c__1);
    trl_g_sum(info->mpicom, 1, &(wrk[0].r), &(wrk[1].r));

    if (!(wrk[0].r >= zero) || !(wrk[0].r <= DBL_MAX)) {
	info->stat = -102;
	return;
    }
    beta[jnd - 1] = sqrt(wrk[0].r);
    tmp1 = alpha[jnd - 1] * alpha[jnd - 1];
    if (jm1 > kept) {
	tmp1 += (beta[jm1 - 1] * beta[jm1 - 1]);
#if TUNED < 2
	(info->flop) += (2 * nrow + 4);
#endif
    } else if (kept > 0) {
	tmp1 += trl_ddot(jm1, beta, c__1, beta, c__1);
#if TUNED < 2
	(info->flop) += (2 * (nrow + kept + 2));
#endif
    }

    if (jm1 == kept) {
	usecgs = 1;
    } else if (jnd >= info->ntot) {
	usecgs = 0;
    } else if (DBL_EPSILON * wrk[0].r >= tmp1) {
	double anorm = 0.0;
	for (i = 0; i < jnd; ++i) {
	    d__1 = fabs(alpha[i]) + fabs(beta[i]);
	    if (d__1 > anorm)
		anorm = d__1;
	}
	usecgs = (beta[jm1] < DBL_EPSILON * anorm * info->ntot);
    } else {
	usecgs = 1;
    }

// whether to perform full re-orthogonalization or extended local
// re-orthogonalization
    if (usecgs != 0) {
	// perform global re-orthogonalization
	nr = info->nrand;
	no = info->north;
	info->stat = ztrl_cgs(info, nrow, v1, ld1, m1, v2, ld2, m2, rr,
			       &beta[jnd - 1], &alpha[jnd - 1],
			       &(info->north), wrk);

#if TUNED < 2
	(info->flop) +=
	    (4 * nrow * ((info->north - no) * jnd + info->nrand - nr) +
	     nrow);
#endif

    } else if (jnd > 1) {
	// perform local re-orthogonalization against two previous vectors
	//printf( " local re-orthogonalization against two previous vectors.\r\n" );
	if (m2 > 1) {
	    qa = &v2[(m2 - 1) * ld2];
	    qb = &v2[(m2 - 2) * ld2];
	} else if (m2 == 1) {
	    qa = v2;
	    qb = &v1[(m1 - 1) * ld1];
	} else {
	    qa = &v1[(m1 - 1) * ld1];
	    qb = &v1[(jm1 - 1) * ld1];
	}
	wrk[0].r = zero;
	wrk[0].i = zero;
	wrk[1].r = zero;
	wrk[1].i = zero;
	for (i = 0; i < nrow; i++) {
	    wrk[0].r += (qa[i].r * rr[i].r) + (qa[i].i * rr[i].i);
	    wrk[0].i += (qa[i].r * rr[i].i) - (qa[i].i * rr[i].r);

	    wrk[1].r += (qb[i].r * rr[i].r) + (qb[i].i * rr[i].i);
	    wrk[1].i += (qb[i].r * rr[i].i) - (qb[i].i * rr[i].r);
	}
	ztrl_g_sum(info->mpicom, 2, &wrk[0], &wrk[2]);
//   ** updating alpha **
//     alpha[jnd-1] += wrk[0];
	alpha[jnd - 1] += wrk[0].r;
	z__1.r = -wrk[0].r;
	z__1.i = -wrk[0].i;
	trl_zaxpy(nrow, z__1, qa, c__1, rr, c__1);
	z__1.r = -wrk[1].r;
	z__1.i = -wrk[1].i;
	trl_zaxpy(nrow, z__1, qb, c__1, rr, c__1);
//     zdotc_( &z__1,&nrow,rr,&c__1,rr,&c__1 );
//     beta[jnd-1]=sqrt(z__1.r);
	d__1 = one / beta[jnd - 1];
	trl_zdscal(nrow, d__1, rr, c__1);
	//if (info->verbose > -1) {
	//  ztrl_g_dot_(info->mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk);
	//  printf("Orthogonality level of v(%d) (after local reothogonalization):\n", jnd+1);
	//}

#if TUNED < 2
	info->flop = info->flop + 9 * nrow;
#endif

    } else {
	// perform local re-orthogonalization against the only vector
	if (m1 == 1) {
	    qa = v1;
	} else {
	    qa = v2;
	}
	trl_zdotc(wrk, nrow, qa, c__1, rr, c__1);
	ztrl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
//     alpha[jnd-1] += wrk[0]..;
	alpha[jnd - 1] += wrk[0].r;

	z__1.r = -wrk[0].r;
	z__1.i = -wrk[0].i;
	trl_zaxpy(nrow, z__1, qa, c__1, rr, c__1);
	d__1 = one / beta[jnd - 1];
	trl_zdscal(nrow, d__1, rr, c__1);

#if TUNED < 2
	info->flop = info->flop + 5 * nrow;
#endif
    }
    // when beta(jnd) is exceedingly small, it should be treated as zero

    if (info->stat == 0) {
	if (beta[jnd - 1] <= DBL_EPSILON * fabs(alpha[jnd - 1])) {
	    beta[jnd - 1] = zero;
	} else if (jnd >= info->ntot) {
	    beta[jnd - 1] = zero;
	}
    }
//
// .. end of trl_orth ..
//
}

////
int ztrl_cgs(trl_info * info, int nrow, trl_dcomplex * v1, int ld1,
	      int m1, trl_dcomplex * v2, int ld2, int m2,
	      trl_dcomplex * rr, double *rnrm, double *alpha, int *north,
	      trl_dcomplex * wrk)
{
//
// Purpose
// =======
// Perform full Gram-Schmidt routine -- orthogonalize a new vector against all existing vectors.
//
// Arguments
// =========
// info   (input) Pointer to structure trl_info_
//         On entry, points to the current TRL_INFO.
//
// nrow   (input) Integer
//         On entry, specifies the number of rows in eigenvectors.
//
// v1     (input) Complex array (ld1,m1)
//         On entry, contains the first part of Lanczos basis computed.
//
// ld1    (input) Integer
//         On entry, specifies the leading dimention of v1.
//
// m1     (input) Integer
//         On entry, specifies the number of Lanczos basis in v1.
//
// v2     (input) Complex array (ld2,m2)
//         On entry, contains the second part of Lanczos basis computed.
//
// ld2    (input) Integer
//         On entry, specifies the leading dimention of v2.
//
// m2     (input) Integer
//         On entry, specifies the number of Lanczos basis in v2.
//
// rr     (input/output) Complex array (nrow)
//         On entry, contains the new Lanczos basis computed.
//         On exit, contains the next Lanczos basis computed after the orthogonalization.
//
// rnrm   (output) double precision
//         On entry, specifies the norm of the current Lanczos basis.
//         On exit, specifies the norm of the new Lanczos basis.
//
// alpha  (input/output) double precision array (m1+m2)
//         On entry, contains the alpha values, on exit, they are updated.
//
// north  (output)
//         On exit, specifies the number of times the full-orthogonalization is applied.
//
// wrk    (workspace) complex array (m1+m2)
//
// ..
// .. local parameters ..
    long c__1 = 1;
    char notrans = 'N';
    double one = 1.0, zero = 0.0;
    const int maxorth = 3;
//
// ..
// .. local variables ..
    trl_dcomplex z__1, z__2;
    int i, j, k, mpicom, myid, nold, irnd, cnt, ierr = 0;
    double tmp, old_rnrm, d__1;
    trl_dcomplex tmp2;
//
// ..
// .. executable statements ..
    mpicom = info->mpicom;
    myid = info->my_pe;
    nold = m1 + m2;
    if (ld1 < nrow || (ld2 < nrow && m2 > 0)) {
	return -201;
    }
    irnd = 0;
    ierr = 0;
    //printf( "nold = %d rnrm = %e\n",nold,*rnrm );
    if (nold > 0) {
	cnt = 0;
	while (cnt <= maxorth) {
	    // compute [v1 v2]'*rr=wrk
	    ztrl_g_dot_(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk);
	    if (m1 > 1) {
		z__1.r = -one;
		z__1.i = zero;
		z__2.r = one;
		z__2.i = zero;
		trl_zgemv(&notrans, nrow, m1, z__1, v1, ld1, wrk, c__1,
			  z__2, rr, c__1);
	    } else if (m1 == 1) {
		z__1.r = -wrk[0].r;
		z__1.i = -wrk[0].i;
		trl_zaxpy(nrow, z__1, v1, c__1, rr, c__1);
	    }
	    if (m2 > 1) {
		z__1.r = -one;
		z__1.i = zero;
		z__2.r = one;
		z__2.i = zero;
		trl_zgemv(&notrans, nrow, m2, z__1, v2, ld2, &wrk[m1],
			  c__1, z__2, rr, c__1);
	    } else if (m2 == 1) {
		z__1.r = -wrk[nold - 1].r;
		z__1.i = -wrk[nold - 1].i;
		trl_zaxpy(nrow, z__1, v2, c__1, rr, c__1);
	    }
	    /* when failing this test, most likely, we should perturb the vector rr.
	       if( irnd == 0 ) {
	       tmp = (ld1*nold)*DBL_EPSILON*max(fabs(*alpha), *rnrm);
	       d__1 = sqrt(wrk[nold-1].r*wrk[nold-1].r + wrk[nold-1].i*wrk[nold-1].i);
	       printf( "sanity check against v[%d]: |v[%d]'*r| = %g <= %g\n", nold, nold, d__1, tmp);
	       if (d__1 > tmp && tmp > zero) {
	       return -202;
	       }
	       // *alpha = *alpha + wrk[nold-1];
	       }
	     */
	    (*north)++;
	    cnt++;
	    trl_zdotc(&tmp2, nold, wrk, c__1, wrk, c__1);
	    tmp = tmp2.r;
	    trl_zdotc(wrk, nrow, rr, c__1, rr, c__1);
	    trl_g_sum(mpicom, 1, &(wrk[0].r), &(wrk[1].r));
	    *rnrm = sqrt(wrk[0].r);

	    old_rnrm = sqrt(wrk[0].r + tmp2.r);
	    //
	    // decisions about whether to re-orthogonalize is based on
	    // relative size of tmp and wrk(1) (R norm square)

	    if (DBL_EPSILON * wrk[0].r > tmp2.r) {
		//if (cnt > 1) {
		//  ztrl_g_dot_(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk);
		//  printf("Orthogonality level of v(%d) (as rr):\n", nold+1);
		//}
		// no need for more orthogonalization
		cnt = maxorth + 1;
	    } else
		if (((cnt > 1
		      && !(tmp2.r <=
			   info->ntot * DBL_EPSILON * DBL_EPSILON *
			   (wrk[0].r + tmp2.r))) ||
		     /*(cnt > 1 && ! (tmp2.r <= info->ntot * DBL_EPSILON * (wrk[0].r+tmp2.r))) || */
		     !(wrk[0].r > DBL_MIN)) && irnd < maxorth) {
		// the input vector is so nearly linear dependent on the
		// existing vectors, we need to perturb it in order to
		// generate a new vector that is orthogonal to the existing
		// ones
		// the perturbation is done in two stages:
		// -- perturbing only one number
		// -- call random_number to generate a whole random vector
		cnt = 0;
		irnd++;
		info->nrand++;
    GetRNGstate();
		if (irnd == 1 && *rnrm > 0.0
		    && (*rnrm) > DBL_EPSILON * old_rnrm) {
		    // old version:  modify only one element
		    // new version:  modify (nrow * epsilon * rnrm / old_rnrm ) elements
		    tmp = unif_rand();
		    i = (int) (nrow * tmp);
		    k = i +
          (int) (fmax2(1.0,
                       (nrow * sqrt(DBL_EPSILON * old_rnrm / *rnrm))));
		    for (j = i; j < k; j++) {
          tmp = unif_rand();
          while (fabs(tmp - 0.5) <= DBL_EPSILON) {
            tmp = unif_rand();
          }
          rr[j].r += (*rnrm) * (tmp - 0.5);
          rr[j].i += (*rnrm) * (tmp - 0.5);
			//rr[j].r += (tmp - 0.5);
			//rr[j].i += (tmp - 0.5);
		    }
		} else {
		    // fill with random numbers produced by intrinsic function
      for (i = 0; i <= myid; i++) {
        tmp = unif_rand();
      }
      i = (int) (nrow * tmp);
      tmp = unif_rand();
		    j = (int) (nrow * tmp);
		    if (i < j) {
			for (k = i; k <= j; k++) {
			    rr[k].r = unif_rand();
			    rr[k].i = unif_rand();
			}
		    } else if (i > j) {
          for (k = j; k <= i; k++) {
            rr[k].r = unif_rand();
            rr[k].i = unif_rand();
          }
		    } else {
          for (k = 0; k < nrow; k++) {
            rr[k].r = unif_rand();
            rr[k].i = unif_rand();
          }
		    }
		}
    PutRNGstate();
		//rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1)
		trl_smooth_zz(nrow, rr);
    }
	}
	// failed to reduce the dot-products between the new vector
	// and the old vectors to small number -- usually an indication of
	// problem with orthogonality of the old vectors
	if (!(wrk[0].r >= tmp))
	    ierr = - 203;
    }
    //
    // normalization
    //
    if (ierr == 0) {
	if (*rnrm > DBL_MIN) {
	    d__1 = one / (*rnrm);
	    trl_zdscal(nrow, d__1, rr, c__1);
	} else {
	    return -204;
	}
    }
    if (irnd > 0)
	*rnrm = zero;
    return ierr;
//
// .. end of trl_cgs ..
//
}
