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
  Xtilde <- out$Xtilde
  Ytilde <- out$Ytilde
  beta_start <- rep(0.5, p)
  lambda <- 2
  f_obj <- crossprod(Ytilde - Xtilde %*% beta_start) / (2 * n) + lambda * sum(abs(beta_start))
  testthat::expect_equal(lasso(Xtilde, Ytilde, beta_start, lambda), f_obj)
  
  # Tests for fitLASSOstandardized()
  testthat::expect_no_error(fitLASSOstandardized(Xtilde, Ytilde, lambda))
  testthat::expect_error(
    fitLASSOstandardized(Xtilde, rep(1, n + 1), lambda),
    "The number of rows in Xtilde should be equal to the length of Ytilde."
  )
  testthat::expect_error(fitLASSOstandardized(Xtilde, Ytilde, -lambda),
                         "lambda should be non-negative.")
  testthat::expect_no_error(fitLASSOstandardized(Xtilde, Ytilde, lambda, beta_start))
  testthat::expect_error(
    fitLASSOstandardized(Xtilde, Ytilde, lambda, rep(0, p + 1)),
    "The length of beta_start should be equal to the number of columns in Xtilde."
  )
  
  # Tests for fitLASSOstandardized_seq()
  lambda_seq <- c(-3, 0, 2, 5, 1)
  testthat::expect_error(
    fitLASSOstandardized_seq(Xtilde, rep(1, n + 1)),
    "The number of rows in Xtilde should be equal to the length of Ytilde."
  )
  testthat::expect_warning(
    fitLASSOstandardized_seq(Xtilde, Ytilde, c(-1, -2)),
    "All values for lambda are less than zero."
  )
  testthat::expect_no_error(out <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq))
  testthat::expect_equal(out$lambda_seq, c(5, 2, 1, 0))
})
