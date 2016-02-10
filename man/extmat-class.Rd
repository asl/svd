\name{extmat-class}
\Rdversion{1.1}
\docType{class}
\alias{extmat-class}
\alias{\%*\%,extmat,extmat-method}
\alias{\%*\%,extmat,matrix-method}
\alias{\%*\%,extmat,numeric-method}
\alias{\%*\%,matrix,extmat-method}
\alias{\%*\%,numeric,extmat-method}
\alias{as.array,extmat-method}
\alias{as.integer,extmat-method}
\alias{as.logical,extmat-method}
\alias{as.matrix,extmat-method}
\alias{as.numeric,extmat-method}
\alias{as.vector,extmat,missing-method}
\alias{crossprod,ANY,extmat-method}
\alias{crossprod,extmat,ANY-method}
\alias{crossprod,extmat,extmat-method}
\alias{crossprod,extmat,missing-method}
\alias{dim,extmat-method}
\alias{length,extmat-method}
\alias{t,extmat-method}
\alias{tcrossprod,ANY,extmat-method}
\alias{tcrossprod,extmat,ANY-method}
\alias{tcrossprod,extmat,extmat-method}
\alias{tcrossprod,extmat,missing-method}

\title{Class \code{"extmat"}}
\description{
  'extmat' is a convenient wrapper which allows one provide just the
  routines which will multiply with matrix and the transposed one
  (e.g. if the matrix is sparse or structured) and allow to use the SVD
  routines of the package. This S4 wrapper allows the use of usual
  matrix operations on such objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{extmat(mul, tmul, nrow, ncol, env = parent.frame())}.
}

\seealso{
  \code{\link{extmat}}
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
\keyword{classes}
