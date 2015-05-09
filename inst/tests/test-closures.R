library(testthat)
library(svd)

as.extmat <- function(m) {
  m <- as.matrix(m)

  extmat(mul = function(x) m %*% x,
         tmul = function(x) x %*% m,
         env = environment(),
         ncol = ncol(m), nrow = nrow(m))
}


t.extmat <- function(emat) {
  emat.ptr <- emat@.xData

  extmat(mul = function(v) .ematmul(emat.ptr, v, transposed = TRUE),
         tmul = function(v) .ematmul(emat.ptr, v, transposed = FALSE),
         nrow = ncol(emat), ncol = nrow(emat))
}

test_that("Transposition (using closures) properly work", {
  set.seed(1)
  m <- matrix(rnorm(12), 3, 4)

  e <- as.extmat(m)

  expect_equal(as.matrix(t(e)), t(m))
  expect_equal(as.matrix(t(t(e))), m)
})
