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
#include <time.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#include "ztrlan.h"
#include "trlan.h"
#include "trlan_i.h"
#include "trlcore_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"
#include "ztrlcore_i.h"
#include "ztrl_comm_i.h"
#include "ztrlan_i.h"

////
void
ztrlan(ztrl_matprod op, trl_info * info, int nrow, int mev, double *eval,
       trl_dcomplex * evec, int lde, trl_dcomplex * misc, int nmis,
       double *dwrk, int ldwrk, void *lparam)
{
//
// Purpose: Top (user) level routines
// ========
// A thick-restart Lanczos routine for computing eigenvalues and
// eigenvectors of a complex hermisian operator/matrix (A).
// -- only accept one input vector, the input vector is expected
//    to be stored in the (nec+1)st column of EVEC.
// -- it extends the Lanczos basis one vector at a time.
// -- orthogonality among the Lanczos vectors are maintained using
//    full re-orthogonalization when necessary.
//
// Requirements:
// =============
// 1) User supplies OP with the specified interface.
// 2) If (info->nec>0), evec(1:nrow, 1:info->nec) must contain valid
//    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
//    These eigenpairs are assumed to have zero residual norms inside
//    TRLAN.
// 3) lde >= nrow.
// 4) The arrays evec and eval must be large enough to store the
//    solutions, i.e., mev >= info->ned and mev >= info->nec.
// 5) The array wrk may be of arbitrary size.  Internally, the workspace
//    size is
//        nrow*max(0,info->ned-size(evec,2))+maxlan*(maxlan+10)
//
//    If wrk is smaller than this, trlan routine will allocate additional
//    workspace to accommodate.
//
// Arguments:
// ==========
// op      (input) function pointer
//         On entry, points to a function that comptues op(X) == A*X,
//         when given a set of vectors X.
//         The operator that defines the eigenvalue problem is expected to
//         have the following interface
//          void op(nrow, ncol, xin, ldx, yout, ldy)
//            nrow  (input) integer
//                   On entry, specifies the number of rows in xin and yout.
//            ncol  (input) integer
//                   On entry, specifies, the number of columns in Xin and
//                   Yout.
//            xin   (input) double precision vector of length (ldx*ncol)
//                   On entry, contatins the input vector to be multiplied.
//                   It consists of ncol column vectors with each column
//                   stored in consecutive order.
//            ldx   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   Xin, i.e., the i-th column vector starts with element
//                   (i-1)*ldx+1 and ends with element (i-1)*ldx+nrow in Xin.
//            yout  (output) double precision vector of length (ldy*ncol)
//                   On exit, contains the result array, i.e., it stores the
//                   result of matrix-vector multiplications.
//            ldy   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   yout, i.e., the i-th column vector starts with element
//                   (i-1)*ldy+1 in Yout and ends with element (i-1)*ldy+nrow.
//
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN
//
// nrow    (input) integer
//          On entry, specifies the number of rows that is on this processor.
//
// mev     (input) integer
//          On entry, specifies the number of Ritz pairs, that can be stored in
//          eval and evec.
//
// eval    (output) double precision vector of length (mev)
//          On exist, stores the eigenvalues
//
// evec    (output) double precision vector of lenvth (nrow*mev)
//          On exit, stores the eigenvectors
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//
// misc     (workspace) complex double vector of length (nmis)
//          If it is provided and there is enough space, the residual norm of
//          the converged eigenpairs will be stored at wrk(1:info->nec) on exit.
//          Internally, it is used as the workspace (size of kept<maxlan) for
//          orthogonalization. It is also used to multiply the ritz vector with
//          lanczos basis (at most nrow*maxlan).
//
// lwrk    (optional) integer
//          On entry, specifies, the size of misc.
//
// dwrk    (workspace) double precision vector of length (ldwrk)
//          On entry, provides a workspace (at leaset 5*jnd, where jnd<maxlan)
//          to compute eigenvalues and eigenvectors of tri-diagonal matrix. It is
//          also used to apply the rotation matrix to the eigenvector computed
//          (at most maxlan*maxlan).
//
// ldwrk   (input) integer
//          On entry, specified the size of workspace.
//
//
// ..
// .. local scalars ..
  clock_t clk1;
  int ii, nbas, ldb;
//
// ..
// .. executables statements ..
  clk1 = clock();
  info->clk_tot = clk1;
  if (info->ned > mev) {
    Rprintf("info->ned (%d) is larger than mev (%d) reducing info->ned "
            "to %d\n", info->ned, mev, mev);
    info->ned = mev;
  }
  //
  // there is nothing to do if there is no more eigenvalue to compute
  if (info->ned <= info->nec || info->ned <= 0)
    goto end;
  //
  info->stat = 0;
  trl_clear_counter(info, nrow, mev, lde);
  if (info->stat != 0)
    goto end;
  //
  // Internally, the workspace is broken into two parts
  // one to store (maxlan+1) Lanczos vectors, and the other to
  // store all others (size maxlan*(maxlan+ned+10))
  // The next If-block decides how the two arrays are mapped.
  //
  //printf( "nmis=%d\n",nmis );
  memset(misc, 0, nmis * sizeof(trl_dcomplex));
  //
  // make sure every process is successful so far
  //printf( "synchronizing.\n" );
  ii = trl_sync_flag(info->mpicom, info->stat);
  info->stat = ii;
  if (ii != 0)
    goto end;
  //
  // open log and checkpoint files
  trl_open_logfile(info);
  if (info->verbose > 0) {
    nbas = imax2(1, info->maxlan - mev + 1);
    ldb = ((nrow + 3) / 4) * 4;
    if ((ldb % 4096) == 0)
      ldb = ldb + 4;
    trl_time_stamp(info->log_fp);
    trl_print_setup(info, nbas * ldb, nmis, ldwrk);
  }
  //
  // call trlanczos to do the real work
  ztrlanczos(op, info, nrow, mev, eval, evec, lde, misc, nmis, dwrk, ldwrk, lparam);
  //
  // close log and checkpoint files
  trl_close_logfile(info);
  //
  // DONE, reclaim the space allocated
end:
  clk1 = clock();
  if (clk1 < info->clk_tot) {
    info->tick_t +=
      (info->clk_max -
       info->clk_tot) / (double) (info->clk_rate);
    info->tick_t +=
      (info->clk_max + clk1) / (double) (info->clk_rate);
    info->clk_tot = 0;
  } else if (info->clk_tot < 0 && clk1 >= 0) {
    info->tick_t -= info->clk_tot / (double) (info->clk_rate);
    info->tick_t += clk1 / (double) (info->clk_rate);
    info->clk_tot = 0;
  } else {
    info->tick_t  += (clk1 - info->clk_tot) / (double) (info->clk_rate);
    info->clk_tot  = 0;
  }
  /*
    ii = clock();
    if (ii >= info->clk_tot) {
    info->clk_tot = ii - info->clk_tot;
    } else if (ii < 0) {
    //assuming ( info->clk_tot > 0 ), otherwise ?? wrap-around twice ??
    info->tick_t +=
    (info->clk_max - info->clk_tot) / (double) (info->clk_rate);
    info->clk_tot = info->clk_max + ii;
    } else {
    //assuming ( info->clk_tot < 0 ),
    info->tick_t +=
    (info->clk_max + info->clk_tot) / (double) (info->clk_rate);
    info->clk_tot = ii;
    }
  */
  return;
//
// .. end of trlan_ ..
//
}

////
void
ztrl_check_ritz(ztrl_matprod op, trl_info * info, int nrow, int ncol,
                trl_dcomplex * rvec, int ldrvec, double *alpha,
                int *check, double *beta, double *eval,
                trl_dcomplex * wrk, int lwrk, void *lparam)
{
//
// Purpose:
// ========
// Performs a standard check on the computed Ritz pairs.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// nrow     (input) integer
//           On entry, specifies the problem size.
//
// ncol     (input) integer
//           On entry, specifies the number of Ritz values computed.
//
// rvec     (input) double precision array of dimension (nrow,ncol)
//           On entry, specifies the array storing the Ritz vectors.
//
// alpha    (input) double precision array of dimension (ncol)
//           On entry, contains he Ritz values computed.
//
// beta     (optional) double precision array of dimension (ncol)
//           If provided, contaions the residual norms returned from a Lanczos routine.
//
// eval     (optional) double precision array of dimension (ncol)
//           If provided, contains the actual eigenvalues.
//
// lwrk     (optional) integer
//           If provided, specifies the size of workspace provided.
//
// wrk      (optional) double precision array of size(lwrk)
//           If provided, double precidion workspace.
//
//
// ..
// .. local parameters ..
  int i__1 = 1;
  long c__1 = 1;
  trl_dcomplex z__1;
//
// ..
// .. local variables ..
// aq -- store the result of matrix-vector multiplication, size nrow
// rq -- store the Rayleigh-Quotient and the residual norms
// gsumwrk -- workspace left over for trl_g_sum to use dimension of the input arrays
  trl_dcomplex *aq, *gsumwrk;
  double *rq, *res, *err, *dsumwrk;
  FILE *fp;
  int i, aqi, gsumwrki, icheck;
  double gapl, gapr;
//
// ..
// .. executable statements ..
  if (ncol <= 0)
    return;			// there is nothing to do
  // figure out whether it is necessary to allocate new workspace
  *check = 0;
  aqi = 0;
  gsumwrki = 0;
  if (lwrk >= nrow + ncol) {
    aq = &wrk[0];
    gsumwrk = &wrk[nrow];
  } else if (lwrk >= ncol) {
    gsumwrk = &wrk[0];
    aq = Calloc(nrow, trl_dcomplex);
    if (aq == NULL) {
      error("TRL_CHECK_RITZ: Failed to allocated workspace AQ");
    }
    aqi = 1;
  } else if (lwrk >= ncol) {
    gsumwrk = wrk;
    aq = Calloc(nrow, trl_dcomplex);
    if (aq == NULL) {
      error("TRL_CHECK_RITZ: Failed to allocated workspace AQ");
    }
    aqi = 1;
  } else {
    // WRK not provided -- allocate space for AQ and RQ,
    // gsumwrk points to the last third of RQ
    aq = Calloc(nrow, trl_dcomplex);
    if (aq == NULL) {
      error("TRL_CHECK_RITZ: Failed to allocated workspace AQ");
    }
    aqi = 1;
    gsumwrk = Calloc(ncol, trl_dcomplex);
    if (gsumwrk == NULL) {
      error("TRL_CHECK_RITZ: Failed to allocate workspace GSUMWRK.\n");
    }
    gsumwrki = 1;
  }
  dsumwrk = Calloc(ncol, double);
  memset(aq, 0, nrow * sizeof(trl_dcomplex));
  memset(gsumwrk, 0, ncol * sizeof(trl_dcomplex));
  memset(dsumwrk, 0, ncol * sizeof(double));
  //
  // go through each Ritz pair one at a time, compute Rayleigh
  // quotient and the corresponding residual norm
  rq = Calloc(3 * ncol, double);
  res = Calloc(ncol, double);
  err = Calloc((2 * ncol), double);
  for (i = 0; i < ncol; i++) {
    op(&nrow, &i__1, &rvec[i * ldrvec], &ldrvec, aq, &nrow, lparam);
    // Rayleigh quotient -- assuming rvec(:,i) has unit norm
    trl_zdotc(&z__1, nrow, &rvec[i * ldrvec], c__1, aq, c__1);
    //ztrl_g_sum(info->mpicom, 1, &z__1, gsumwrk);
    rq[i] = z__1.r;
    trl_g_sum(info->mpicom, 1, &rq[i], dsumwrk);
    zdaxpy_(nrow, -z__1.r, &rvec[i * ldrvec], aq);
    trl_zdotc(&z__1, nrow, aq, c__1, aq, c__1);
    res[i] = z__1.r;
    trl_zdotc(&z__1, nrow, &rvec[i * ldrvec], c__1, &rvec[i * ldrvec],
              c__1);
  }
  trl_g_sum(info->mpicom, ncol, res, dsumwrk);
  for (i = 0; i < ncol; i++) {
    res[i] = sqrt(res[i]);
  }
  //
  // compute the error estimate based on computed residual norms and the
  // Ritz values
  gapl = alpha[ncol - 1] - alpha[0];
  for (i = 0; i < ncol - 1; i++) {
    gapr = alpha[i + 1] - alpha[i];
    gapl = fmin2(gapl, gapr);
    if (res[i] >= gapl) {
      err[i] = res[i];
    } else {
      err[i] = res[i] * res[i] / gapl;
    }
    gapl = gapr;
  }
  if (res[ncol - 1] >= gapl) {
    err[ncol - 1] = res[ncol - 1];
  } else {
    err[ncol - 1] = res[ncol - 1] * res[ncol - 1] / gapl;
  }
  // if writing to stdout, only PE 0 does it
  if (info->log_fp == NULL) {
    trl_reopen_logfile(info);
  }
  if (info->my_pe <= 0) {
    if (info->stat != 0) {
      *check = -4;
    }
    for (i = 0; i < ncol; i++) {
      if (err[i] < DBL_EPSILON * alpha[ncol - 1]) {
        err[i] = DBL_EPSILON * alpha[ncol - 1];
      }
    }
    // print out the information
    Rprintf("TRL_CHECK_RITZ: \n");
    Rprintf(
            "           Ritz value       res norm   res diff  est error  diff w rq  act. error\n");
    if (beta != NULL && eval != NULL) {
      for (i = 0; i < ncol; i++) {
        icheck = 0;
        Rprintf(
                "%21.14f    %11.3e%11.3e%11.3e%11.3e %11.3e%11.3e\n",
                alpha[i], res[i], beta[i] - res[i], err[i],
                rq[i] - alpha[i], eval[i] - alpha[i], eval[i]);
        //
        // check the accuracy of results..
        if (fabs(beta[i] - res[i]) > 0.00001) {
          Rprintf(
                  " res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
                  i, beta[i], res[i], fabs(beta[i] - res[i]));
          *check = *check - 1;
          icheck++;
        }
        if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
          Rprintf(
                  " diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
                  i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
                  nrow * nrow * info->tol);
          *check = *check - 1;
          icheck++;
        }
        if (fabs(eval[i] - alpha[i]) > 10 * nrow * nrow * info->tol ||
            fabs(eval[i] - alpha[i]) > 10 * err[i]) {
          Rprintf(
                  " act. error[%d] = |exact-alpha| = %5.3e - %5.3e = %5.3e > 10*nrow*nrow*tau =%5.3e or 10*est err = %5.3e\n",
                  i, eval[i], alpha[i], fabs(eval[i] - alpha[i]),
                  10 * nrow * nrow * info->tol, 10 * err[i]);
          *check = *check - 1;
          icheck++;
        }
        if (icheck > 0) {
          Rprintf("\n");
        }
      }

    } else if (beta != NULL) {
      for (i = 0; i < ncol; i++) {
        Rprintf("%21.14f    %11.3e%11.3e%11.3e%11.3e\n",
                alpha[i], res[i], beta[i] - res[i], err[i],
                rq[i] - alpha[i]);
        //
        // check the accuracy of results..
        if (fabs(beta[i] - res[i]) > 0.00001) {
          Rprintf(
                  " res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
                  i, beta[i], res[i], fabs(beta[i] - res[i]));
          *check = *check - 1;
          icheck++;
        }
        if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
          Rprintf(
                  " diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
                  i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
                  nrow * nrow * info->tol);
          *check = *check - 1;
          icheck++;
        }
        if (icheck > 0) {
          Rprintf("\n");
        }
      }

    } else if (eval != NULL) {
      for (i = 0; i < ncol; i++) {
        Rprintf(
                "%21.14f     %11.3e           %11.3e%11.3e%11.3e%11.3e\n",
                alpha[i], res[i], err[i], rq[i] - alpha[i],
                eval[i] - alpha[i], eval[i]);
      }
    } else {
      for (i = 0; i < ncol; i++) {
        Rprintf("%21.14f    %11.5e           %11.3e%11.3e\n",
                alpha[i], res[i], err[i], rq[i] - alpha[i]);
      }
    }
  }
  if (info->nec < info->ned)
    *check = 1;
  Free(res);
  Free(err);
  Free(rq);
  Free(dsumwrk);
  if (aqi > 0) {
    Free(aq);
  }
  if (gsumwrki > 0) {
    Free(gsumwrk);
  }
  trl_close_logfile(info);
//
// .. end of trl_check_ritz_
//
}
