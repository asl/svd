\name{extmat}
\alias{is.extmat}
\alias{extmat.ncol}
\alias{extmat.nrow}

\title{External matrices operations.}

\description{
  A set of routines to operate on "external" matrices.
}

\usage{
is.extmat(X)
extmat.ncol(X)
extmat.nrow(X)
}

\arguments{
  \item{X}{matrix to operate on.}
}

\details{
  These routines checks whether the given external pointer actually
  points to "external matrix" structure and allow to extract the number
  of columns and rows respectively.
}
