library(testthat)
library(svd)

as.extmat <- function(m) {
  m <- as.matrix(m)

  extmat(mul = function(x) m %*% x,
         tmul = function(x) x %*% m,
         env = environment(),
         ncol = ncol(m), nrow = nrow(m))
}


test_that("Multiplication on vector works ok", {
  set.seed(1)
  m <- matrix(rnorm(12), 3, 4)
  m2 <- matrix(rnorm(12), 3, 4)

  e <- as.extmat(m)

  v <- rnorm(ncol(e))
  expect_equal(e %*% v, m %*% v)
  expect_equal(tcrossprod(e, v), m %*% v)

  w <- rnorm(nrow(e))

  expect_equal(w %*% e, w %*% m)
  expect_equal(crossprod(e, w), t(m) %*% w)

  expect_equal(crossprod(e), crossprod(m))
  expect_equal(tcrossprod(e), tcrossprod(m))

  expect_equal(tcrossprod(e, m2), tcrossprod(m, m2))
  expect_equal(crossprod(e, m2), crossprod(m, m2))
  # expect_equal(crossprod(e, as.extmat(m2)), crossprod(m, m2))
  # expect_equal(tcrossprod(e, as.extmat(m2)), tcrossprod(m, m2))


  expect_equal(e %*% t(m2), m %*% t(m2))
  expect_equal(t(m2) %*% e, t(m2) %*% e)

  expect_equal(e %*% as.extmat(t(m2)), m %*% t(m2))
})

test_that("Empty result", {
  set.seed(1)
  m <- matrix(rnorm(12), 3, 4)
  m2 <- matrix(rnorm(12), 3, 4)

  e <- as.extmat(m)
  expect_equal(e %*% matrix(0, ncol(e), 0), matrix(0, nrow(e), 0))
  expect_equal(matrix(0, 0, nrow(e)) %*% e, matrix(0, 0, ncol(e)))
})
