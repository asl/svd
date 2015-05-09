#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   'svd' R package
#   Copyright (c) 2015 Anton Korobeynikov <asl@math.spbu.ru>
#
#   This program is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#   PURPOSE.  See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program; if not, write to the
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
#   MA 02139, USA.

#   Routines for external matrix stuff

extmat.ncol <- function(X)
  .extmat.ncol(X@.xData)

extmat.nrow <- function(X)
  .extmat.nrow(X@.xData)

.extmat.ncol <- function(X) {
  .Call("extmat_ncol", X)
}

.extmat.nrow <- function(X) {
  .Call("extmat_nrow", X)
}

is.extmat <- function(X)
  is(X, "extmat") && .is.extmat(X@.xData)

.is.extmat <- function(X) {
  .Call("is_extmat", X)
}

ematmul <- function(emat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("ematmul", emat@.xData, v, transposed);
}

.ematmul <- function(emat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("ematmul_unchecked", emat, v, transposed);
}

extmat <- function(mul, tmul, nrow, ncol,
                   env = parent.env()) {
  new("extmat",
      .Call("initialize_rextmat",
            match.fun(mul), match.fun(tmul),
            as.integer(nrow), as.integer(ncol),
            env))
}

# S4 weirdness
setClass("extmat", contains = "externalptr")
  
setMethod("as.matrix", signature(x = "extmat"), function(x) as(x, "matrix"))
setMethod("as.array",  signature(x = "extmat"), function(x) as(x, "matrix"))
as.array.extmat <- as.matrix.extmat <- function(x, ...) as(x, "matrix")

setMethod("as.vector", signature(x = "extmat", mode = "missing"),
          function(x, mode) as.vector(as(x, "matrix"), mode))
as.vector.extmat <- function(x, mode) as.vector(as(x, "matrix"), mode)

setMethod("as.numeric", signature(x = "extmat"),
          function(x, ...) as.numeric(as.vector(x)))
setMethod("as.integer", signature(x = "extmat"),
          function(x, ...) as.integer(as.vector(x)))
setMethod("as.logical", signature(x = "extmat"),
          function(x, ...) as.logical(as.vector(x)))

setAs("extmat", "matrix", function(from) from %*% diag(nrow = ncol(from)))

setMethod("dim", signature(x = "extmat"),
          function(x) c(.extmat.nrow(x@.xData), .extmat.ncol(x@.xData)), valueClass = "integer")
setMethod("length", "extmat", function(x) prod(dim(x)))

setMethod("%*%", signature(x = "extmat", y = "numeric"),
          function(x, y) {
            dim(y) <-
              if (ncol(x) == (n <- length(y))) c(n, 1L) else c(1L, n)
            x %*% y
          })
setMethod("%*%", signature(x = "numeric", y = "extmat"),
          function(x, y) {
            dim(x) <-
              if (nrow(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
            x %*% y
          })

setMethod("%*%", signature(x = "extmat", y = "matrix"),
          function(x, y) {
            if (nrow(y) != ncol(x))
              stop("non-conformable arguments")
            apply(y, 2, .ematmul, emat = x@.xData, transposed = FALSE)
          })
setMethod("%*%", signature(x = "matrix", y = "extmat"),
          function(x, y) {
            if (nrow(y) != ncol(x))
              stop("non-conformable arguments")
            t(apply(x, 1, .ematmul, emat = y@.xData, transposed = TRUE))
          })

# t(m) %*% y
setMethod("crossprod", signature(x = "extmat", y = "matrix"),
          function(x, y) {
            if (nrow(y) != nrow(x))
              stop("non-conformable arguments")
            apply(y, 2, .ematmul, emat = x@.xData, transposed = TRUE)
          })
setMethod("crossprod", signature(x = "extmat", y = "numeric"),
          function(x, y) {
            dim(y) <-
              if (nrow(x) == (n <- length(y))) c(n, 1L) else c(1L, n)
            crossprod(x, y)
          })
# t(m) %*% m
setMethod("crossprod", signature(x = "extmat", y = "missing"),
          function(x, y) {
            # FIXME: get rid of as.matrix, calculate column-by-column
            crossprod(x, as.matrix(x))
          })

# t(t(m) %*% x)
setMethod("crossprod", signature(x = "matrix", y = "extmat"),
          function(x, y) {
            if (nrow(y) != nrow(x))
              stop("non-conformable arguments")
            t(apply(x, 2, .ematmul, emat = y@.xData, transposed = TRUE))
          })
setMethod("crossprod", signature(x = "numeric", y = "extmat"),
          function(x, y) {
            dim(x) <-
              if (nrow(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
            crossprod(x, y)
          })

# m %*% t(y)
setMethod("tcrossprod", signature(x = "extmat", y = "matrix"),
          function(x, y) {
            if (ncol(y) != ncol(x))
              stop("non-conformable arguments")
            apply(y, 1, .ematmul, emat = x@.xData, transposed = FALSE)
          })
setMethod("tcrossprod", signature(x = "extmat", y = "numeric"),
          function(x, y) {
            dim(y) <-
              if (nrow(x) == (n <- length(y))) c(n, 1L) else c(1L, n)
            tcrossprod(x, y)
          })
# m %*% t(m)
setMethod("trossprod", signature(x = "extmat", y = "missing"),
          function(x, y) {
            # FIXME: get rid of diag, calculate column-by-column
            tcrossprod(x, as.matrix(x))
          })

# t(m %*% t(x))
setMethod("tcrossprod", signature(x = "matrix", y = "extmat"),
          function(x, y) {
            if (ncol(y) != ncol(x))
              stop("non-conformable arguments")
            t(apply(x, 1, .ematmul, emat = y@.xData, transposed = FALSE))
          })
setMethod("tcrossprod", signature(x = "numeric", y = "extmat"),
          function(x, y) {
            dim(x) <-
              if (ncol(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
            tcrossprod(x, y)
          })
