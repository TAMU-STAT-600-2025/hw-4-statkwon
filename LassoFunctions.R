# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y) {
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  # [ToDo] Center and scale X
  n <- nrow(X)
  p <- ncol(X)
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, n, p, byrow = TRUE)
  weights <- sqrt(colSums(Xcentered^2) / n)
  Xtilde <- Xcentered %*% diag(1 / weights)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(
    Xtilde = Xtilde,
    Ytilde = Ytilde,
    Ymean = Ymean,
    Xmeans = Xmeans,
    weights = weights
  ))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda) {
  return(sign(a) * max(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda) {
  n <- nrow(Xtilde)
  f_obj <- crossprod(Ytilde - Xtilde %*% beta) / (2 * n) + lambda * sum(abs(beta))
  return(f_obj)
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde,
                                 Ytilde,
                                 lambda,
                                 beta_start = NULL,
                                 eps = 0.001) {
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (n != length(Ytilde)) {
    stop("The number of rows in Xtilde should be equal to the length of Ytilde.")
  }
  
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) {
    stop("lambda should be non-negative.")
  }
  
  #[ToDo]  Check for starting point beta_start.
  if (is.null(beta_start)) {
    # If none supplied, initialize with a vector of zeros.
    beta_start <- rep(0, p)
  } else {
    # If supplied, check for compatibility with Xtilde in terms of p
    if (length(beta_start) != p) {
      stop("The length of beta_start should be equal to the number of columns in Xtilde.")
    }
  }
  
  #[ToDo]  Coordinate-descent implementation.
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  beta_old <- beta_start
  beta <- beta_old
  r <- Ytilde - Xtilde %*% beta
  while (TRUE) {
    for (j in 1:p) {
      beta[j] <- soft(beta_old[j] + crossprod(Xtilde[, j], r) / n, lambda)
      r <- r + Xtilde[, j] * (beta_old[j] - beta[j])
    }
    fmin_old <- lasso(Xtilde, Ytilde, beta_old, lambda)
    fmin <- lasso(Xtilde, Ytilde, beta, lambda)
    if (fmin_old - fmin < eps) {
      break
    }
    beta_old = beta
  }
  
  # Return
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

calculate_lambda_max <- function(Xtilde, Ytilde) {
  n <- nrow(Xtilde)
  return(max(crossprod(Xtilde, Ytilde) / n))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde,
                                     Ytilde,
                                     lambda_seq = NULL,
                                     n_lambda = 60,
                                     eps = 0.001) {
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if (n != length(Ytilde)) {
    stop("The number of rows in Xtilde should be equal to the length of Ytilde.")
  }
  
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  if (!is.null(lambda_seq)) {
    # If lambda_seq is supplied, only keep values that are >= 0,
    # and make sure the values are sorted from largest to smallest.
    lambda_seq <- sort(lambda_seq[lambda_seq >= 0], decreasing = TRUE)
    # If none of the supplied values satisfy the requirement,
    # print the warning message and proceed as if the values were not supplied.
    if (length(lambda_seq) == 0) {
      warning("All values for lambda are less than zero.")
      lambda_max <- calculate_lambda_max(Xtilde, Ytilde)
      lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
    }
  } else {
    # If lambda_seq is not supplied, calculate lambda_max
    # (the minimal value of lambda that gives zero solution),
    # and create a sequence of length n_lambda as
    lambda_max <- calculate_lambda_max(Xtilde, Ytilde)
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda
  # (make sure supplied eps is carried over).
  # Use warm starts strategy discussed in class for setting the starting values.
  n_lambda <- length(lambda_seq)
  beta_mat <- matrix(nrow = p, ncol = n_lambda)
  beta_start <- rep(0, p)
  fmin_vec <- vector(mode = "numeric", length = n_lambda)
  
  for (i in 1:n_lambda) {
    out <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_start, eps)
    beta_mat[, i] <- out$beta
    beta_start <- beta_mat[, i]
    fmin_vec[i] <- out$fmin
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(
    lambda_seq = lambda_seq,
    beta_mat = beta_mat,
    fmin_vec = fmin_vec
  ))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,
                     Y,
                     lambda_seq = NULL,
                     n_lambda = 60,
                     eps = 0.001) {
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  out1 <- standardizeXY(X, Y)
  
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  out2 <- fitLASSOstandardized_seq(out1$Xtilde, out1$Ytilde, lambda_seq, n_lambda, eps)
  lambda_seq <- out2$lambda_seq
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  beta_mat <- out1$weights * out2$beta_mat
  beta0_vec <- out1$Ymean - out1$Xmeans %*% beta_mat
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(
    lambda_seq = lambda_seq,
    beta_mat = beta_mat,
    beta0_vec = beta0_vec
  ))
}

# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,
                    Y,
                    lambda_seq = NULL,
                    n_lambda = 60,
                    k = 5,
                    fold_ids = NULL,
                    eps = 0.001) {
  n <- nrow(X)
  
  # [ToDo] Fit Lasso on original data using fitLASSO
  out1 <- fitLASSO(X, Y, lambda_seq, n_lambda, eps)
  lambda_seq <- out1$lambda_seq
  n_lambda <- length(lambda_seq)
  beta_mat <- out1$beta_mat
  beta0_vec <- out1$beta0_vec
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)) {
    fold_ids <- sample(rep(1:k, length.out = n), size = n)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  cv_folds <- matrix(nrow = k, ncol = n_lambda)
  for (fold in 1:k) {
    Xtrain <- X[fold_ids != fold, ]
    Ytrain <- Y[fold_ids != fold]
    
    Xtest <- X[fold_ids == fold, ]
    Ytest <- Y[fold_ids == fold]
    
    out2 <- fitLASSO(Xtrain, Ytrain, lambda_seq, n_lambda, eps)
    
    n_fold <- nrow(Xtrain)
    for (i in 1:n_lambda) {
      cv_folds[fold, i] <- crossprod(Ytest - out2$beta0_vec[i] - Xtest %*% out2$beta_mat[, i]) / n_fold
    }
  }
  
  cvm <- colMeans(cv_folds)
  cvse <- apply(cv_folds, 2, \(cv_j) {
    sd(cv_j) / sqrt(k)
  })
  
  # [ToDo] Find lambda_min
  lambda_min_idx <- which.min(cvm)
  lambda_min <- lambda_seq[lambda_min_idx]
  
  # [ToDo] Find lambda_1SE
  lambda_max_idx <- which.max(ifelse(cvm <= cvm[lambda_min_idx] + cvse[lambda_min_idx], cvm, 0))
  lambda_1se <- lambda_seq[lambda_max_idx]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(
    list(
      lambda_seq = lambda_seq,
      beta_mat = beta_mat,
      beta0_vec = beta0_vec,
      fold_ids = fold_ids,
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      cvm = cvm,
      cvse = cvse
    )
  )
}
