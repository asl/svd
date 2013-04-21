\name{extmat}
\alias{is.extmat}
\alias{extmat.ncol}
\alias{extmat.nrow}
\alias{extmat}
\alias{ematmul}

\title{External matrices operations.}

\description{
  A set of routines to operate on "external" matrices.
}

\usage{
is.extmat(X)
extmat.ncol(X)
extmat.nrow(X)
extmat(mul, tmul, nrow, ncol, env = .GlobalEnv)
ematmul(emat, v, transposed = FALSE)
}

\arguments{
  \item{X, emat}{matrix to operate on}
  \item{mul}{function performing the multiplication of matrix to vector}
  \item{tmul}{function performing the multiplication of transposed matrix to vector}
  \item{nrow}{number of rows of the matrix}
  \item{ncol}{number of columns of the matrix}
  \item{env}{environment, where matrix-vector multiplication function
    call is evaluated in}
  \item{transposed}{logical, if 'TRUE' the multiplication is performed
    with the transposed matrix.}
  \item{v}{vector to multiply with.}
}

\details{
  These routines checks whether the given external pointer actually
  points to "external matrix" structure and allow to extract the number
  of columns and rows respectively.

  'extmat' is a convenient wrapper which allows one provide just the
  routines which will multiply with matrix and the transposed one
  (e.g. if the matrix is sparse or structured) and allow to use the SVD
  routines of the package
}

\examples{
\dontrun{
library(Matrix)
f <- function(v) as.numeric(A \%*\% v) # Convert Matrix object back to vector
tf <- function(v) as.numeric(tA \%*\% v) # Convert Matrix object back to vector

e <- new.env()
assign("A", USCounties, e)
assign("tA", t(USCounties), e)
environment(f) <- e
environment(tf) <- e

m <-extmat(f, tf, nrow(USCounties), ncol(USCounties))
system.time(v1 <- propack.svd(m, neig = 10))
#   user  system elapsed 
#  0.252   0.007   0.259 
system.time(v2 <- propack.svd(as.matrix(USCounties), neig = 10))
#   user  system elapsed 
#  8.563   0.027   8.590 
}
}
