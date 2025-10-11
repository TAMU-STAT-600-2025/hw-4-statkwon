source("LassoFunctions.R")

testthat::test_that("Test for LassoFunctions.R", {
  # Generate mocking data
  set.seed(0)
  n <- 20
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p)
  epsilon <- rnorm(n)
  Y <- X %*% beta + epsilon
  
  # Tests for standardizeXY()
  out <- standardizeXY(X, Y)
  testthat::expect_equal(diag(crossprod(out$Xtilde) / n), rep(1, p))
  testthat::expect_equal(colMeans(out$Xtilde), rep(0, p))
  testthat::expect_equal(mean(out$Ytilde), 0)
  
  # Tests for soft()
  testthat::expect_equal(soft(5, 3), 2)
  testthat::expect_equal(soft(-5, 3), -2)
  
  # Test for lasso()
  beta_start <- rep(0.5, p)
  lambda <- 2
  f_obj <- crossprod(out$Ytilde - out$Xtilde %*% beta_start) / (2 * n) + lambda * sum(abs(beta_start))
  testthat::expect_equal(lasso(out$Xtilde, out$Ytilde, beta_start, lambda), f_obj)
})
