#ifndef __ZTRLAN_H
#define __ZTRLAN_H
#include "trlan.h" /* trl_info */

#include <R_ext/Complex.h>
typedef Rcomplex trl_dcomplex;

/**
   Prototype matrix-vector multiplication function.  Useful if you have a
   matrix-vector multiplication function defined in Fortran.

   @arg pnrow Pointer to integer value nrow, the number of rows locally on
   this processor.

   @arg pncol Pointer to integer value ncol, the number of columns in
   vectors x and y.

   @arg x Pointer to the elements of input vectors x.  It is assumed to
   have [*pldx * *pncol] elements, with the ith column starting at element
   *pldx * i.

   @arg pldx Pointer to integer value ldx, the leading dimension of x.
   Note that *pldx must not be less than *pnrow.

   @arg y Pointer to the elements of output vectors y.  It is assumed to
   have [*pldy * *pncol] elements, with the ith column starting at element
   *pldy * i.

   @arg pldy Pointer to integer value ldy, the leadying dimension of y.
   Note that *pldy must not be less than *pnrow.
*/
typedef void (*ztrl_matprod) (int *pnrow, int *pncol,
                              trl_dcomplex *x, int *pldx,
                              trl_dcomplex *y, int *pldy,
                              void *lparam);

/*
 Top-level routine for solving Hermitian eigenvalue problems.

 A thick-restart Lanczos routine for computing eigenvalues and
 eigenvectors of a complex Hermitian operator/matrix (A).
 -- It only accept one input vector, the input vector is expected
    to be stored in the (nec+1)st column of EVEC.
 -- It extends the Lanczos basis one vector at a time.
 -- Orthogonality among the Lanczos vectors are maintained using
    full re-orthogonalization when necessary.
 -- User must initialize the trl_info object passed to this function.

 Requirements:
 1) User supplies OP with the specified interface.
 2) If (info%nec>0), evec(1:nrow, 1:info%nec) must contain valid
    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
    These eigenpairs are assumed to have zero residual norms inside
    TRLAN.
 3) lde >= nrow.
 4) The arrays evec and eval must be large enough to store the
    solutions, i.e., mev >= info%ned and mev >= info%nec.
 5) The array wrk may be of arbitrary size.  Internally, the workspace
    size is
        nrow*max(0,info%ned-size(evec,2))+maxlan*(maxlan+10)

    If wrk is smaller than this, trlan routine will allocate additional
    workspace to accommodate.

@arg op      (input) function pointer
         It points to a function that comptues op(X) == A*X,
         when given a set of vectors X.

@arg info    (input) [pointer to the structure trl_info]
          It points to the data structure to store the information
          about the eigenvalue problem and the progress of TRLAN

@arg nrow    (input) [integer]
          It specifies the number of rows that is on this processor.

@arg mev     (input) [integer]
          It specifies the number of Ritz pairs, that can be stored in
          eval and evec.

@arg eval    (output) [double precision vector of length (mev)]
          On exist, stores the eigenvalues

@arg evec    (output) [double complex vector of lenvth (lde*mev)]
          On exit, stores the eigenvectors

@arg lde     (input) [integer]
          On entry, specifies the leading dimension of the array evec.  The
          ith vector in evec starts with element lde*i.

@arg misc     (workspace) [complex double vector of length (nmis)]
          If there is enough space, the residual norm of the converged
          eigenpairs will be stored at wrk(1:info%nec) on exit.
          Internally, it is used as the workspace (size of kept<maxlan) for
          orthogonalization.  It is also used to multiply the ritz vector
          with lanczos basis (at most nrow*maxlan).

@arg nmisc    (optional) [integer]
          It specifies the size of misc.

@arg dwrk    (workspace) [double precision vector of length (ldwrk)]
          It provides a workspace (at leaset 5*jnd, where jnd<maxlan)
          to compute eigenvalues and eigenvectors of tri-diagonal matrix. It is
          also used to apply the rotation matrix to the eigenvector computed
          (at most maxlan*maxlan).

@arg ldwrk   (input) [integer]
          It specifies the size of workspace.
*/ 
void
ztrlan(ztrl_matprod op, trl_info * info, int nrow, int mev, double *eval,
       trl_dcomplex * evec, int lde, trl_dcomplex * misc, int nmis,
       double *dwrk, int ldwrk, void *lparam);

/**
 Performs sanity checks on the computed Ritz pairs.

@arg op       (input) [function pointer]
           It points to the matrix-vector multiplication routine.

@arg nrow     (input) [integer]
           It specifies the problem size.

@arg ncol     (input) [integer]
           It specifies the number of Ritz values computed.

@arg rvec     (input) [double complex array of dimension (ldrvec*ncol)]
           It specifies the array storing the Ritz vectors.

@arg alpha    (input) [double precision array of dimension (ncol)]
           On entry, contains he Ritz values computed.

@arg beta     (input) [double precision array of dimension (ncol)]
           It contaions the residual norms returned from a Lanczos routine.

@arg eval     (input) [double precision array of dimension (ncol)]
           It contains the actual eigenvalues.

@arg lwrk     (input) [integer]
           It specifies the size of workspace provided.

@arg wrk      (workspace) [double complex array of size(lwrk)].
*/
void
ztrl_check_ritz(ztrl_matprod op, trl_info * info, int nrow, int ncol,
		trl_dcomplex * rvec, int ldrvec, double *alpha,
		int *check, double *beta, double *eval,
                trl_dcomplex * wrk, int lwrk, void *lparam);
#endif
