#ifndef __ZTRLCORE_H
#define __ZTRLCORE_H
#include "ztrlan.h"
//
////
void zlog_error_state(trl_info * info, int kept, int j1, int j2,
		      int jnd, int nrow, int mev, double *eval,
		      double *alpha, double *alfrot, double *beta,
		      double *betrot, trl_dcomplex * evec,
		      trl_dcomplex * base, trl_dcomplex * qa,
		      trl_dcomplex * qb, trl_dcomplex * rr, char *title,
		      int *iwrk);
//
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
//             On entry, provides a space to store the title of the information 
//             printed out.
//
// iwrk       (workspace) Integer array of dimension ()
//
// 
// 
//
////
void zprint_restart_state(trl_info * info, char *title, int nrow,
			  int mev, double *alpha, double *beta,
			  double *betrot, trl_dcomplex * evec, int lde,
			  double *yy, int kept, int locked, int *iwrk,
			  double *wrk2, int i2, int jml);
//
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
//
////
void zprint_final_state(trl_info * info, char *title, int nrow, int mev,
			double *eval, int lde, double *beta,
			trl_dcomplex * evec, double *yy, int kept,
			int jml);
//
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
// 
//
////
void zwrite_checkpoint(trl_info * info, char *title, int nrow,
		       double *alpha, double *beta, trl_dcomplex * evec,
		       trl_dcomplex * base, int lde, int j1n, int jnd,
		       int ldb, int j2n);
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
//
////
void
ztrlanczos(ztrl_matprod op, trl_info * info, int nrow, int mev, double *eval,
           trl_dcomplex * evec, int lde, trl_dcomplex * wrk, int lwrk,
           double *dwrk, int ldwrk, void *lparam);
// Purpose
// =======
// The actual work routine of restarted Lanczos program for real
// symmetric eigenvalue problems
//
// user may directly invoke this sunroutine but she/he is responsible
// for input correct arguments as needed
//
// 1) info needs to be initialized
// 2) if info%nec>0, the first nec elements of eval and first nec
//    columns of evec must contain valid eigenpairs
// 3) workspace of currect size
//    eval(mev)
//    evec(lde, mev) (lde >= nrow, mev >= ned)
//    base(nrow, info%maxlan-mev+1) (nrow>=nrow, not used if mev>maxlan)
//    wrk(lwrk) minimum amount of memory required by TRLANCZOS is
//    maxlan*(maxlan+10)
// 4) if log files are to be written, the user need to open files on IO
//    unit log_io so that the log gile may be written correctly.
// 5) the user must set the timing variable info%clk_tot and
//    info%clk_max using system_clock function call in order for this
//    subroutine to track timing results correctly
//
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
////
void ztrl_initial_guess(int nrow, trl_dcomplex * evec, int lde, int mev,
			trl_dcomplex * base, int ldb, int nbas,
			double *alpha, double *beta, int *j1, int *j2,
                        trl_info * info, trl_dcomplex * wrk, int lwrk, void *lparam);
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
////
void ztrl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
		       trl_dcomplex * vec1, int ld1, int m1,
		       trl_dcomplex * vec2, int ld2, int m2,
		       trl_dcomplex * wrk, int lwrk, trl_dcomplex * yy2);
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
////
void ztrl_orth(int nrow, trl_dcomplex * v1, int ld1, int m1,
	       trl_dcomplex * v2, int ld2, int m2, trl_dcomplex * rr,
	       int kept, double *alpha, double *beta, trl_dcomplex * wrk,
	       int lwrk, trl_info * info);
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
////
int ztrl_cgs(trl_info * info, int nrow, trl_dcomplex * v1, int ld1,
	     int m1, trl_dcomplex * v2, int ld2, int m2,
	     trl_dcomplex * rr, double *rnrm, double *alpha, int *north,
	     trl_dcomplex * wrk);
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
//
//
// Math kernels for complex-real operations
// 
int zdgemm_(int m, int n, int k, trl_dcomplex * a, int lda, double *b,
	    int ldb, trl_dcomplex * c, int ldc);
int zdgemm2_(int m, int n, int k, trl_dcomplex * a, int lda, double *b,
	     int ldb, trl_dcomplex * c, int ldc);
int zdgemv_(int m, int n, double alpha, trl_dcomplex * a, int lda,
	    double *x, double beta, trl_dcomplex * y);
int zdaxpy_(int n, double a, trl_dcomplex * zx, trl_dcomplex * zy);

#endif
