\name{eigen}
\alias{trlan.eigen}
\alias{ztrlan.eigen}

\title{
  Generic Eigendecomposition of a Matrix
}

\description{
  Compute the set of eigenvalues and eigenvectors decomposition of a
  real or complex rectangular matrix.
}

\usage{
trlan.eigen(X, neig = min(m, n), opts = list(), lambda = NULL, U = NULL)
ztrlan.eigen(X, neig = min(m, n), opts = list(), lambda = NULL, U = NULL)
}

\arguments{
  \item{X}{the matrix to be decomposed. This can be either normal matrix
  or 'external matrix' object (e.g. one, created via 'extmat' function).}
  \item{neig}{number of desired eigentriples}
  \item{opts}{different options for eigensolver. See 'Details' section
  for more information}
  \item{lambda}{set of already computed singular values (used for
  continuation of the decomposition).}
  \item{U}{matrix of already computed eigenvectors (used for
  continuation of the decomposition).}
}

\details{
  These routines provides an interface to state-of-art
  implementation of eigensolver. In particular, nu-TRLAN does the
  thick-restart Lanczos eigendecomposition of a matrix.
  
  'opts' is a list of different options which can be passed to the
  routines. Note that by default more or less suitable values for these
  options are set by the routines automatically.

  The options for nu-TRLAN are:
  \describe{
    \item{kmax}{integer, maximum number of iterations.}
    \item{maxiter}{integer. maximum number of matrix-vector products.}
    \item{tol}{numeric, tolerance level.}
    \item{verbose}{integer, verboseness level.}
  }
}

\value{
  The returned value is a list with components
  \item{d}{a vector containing the eigenvalues of 'X'}
  \item{u}{a matrix whose columns contain the eigenvectors of
    'X'}
}

\references{
  Wu, K. and Simon, H. (2000). \emph{Thick-restart Lanczos method for
  large symmetric eigenvalue problems.} SIAM J. Matrix Anal. Appl. 22, 2, 602-616.

  Yamazaki, I., Bai, Z., Simon, H., Wang, L.-W., and Wu,
  K. (2008). \emph{Adaptive projection subspace dimension for the thick
  restart Lanczos method.} Tech. rep., Lawrence Berkeley National
  Laboratory, University of California, One Cyclotron road, Berkeley, California 94720.  
  
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
  SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\keyword{array}
\keyword{algebra}
