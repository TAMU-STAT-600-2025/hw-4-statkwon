# Load the riboflavin data
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
? riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
out <- fitLASSO(X, Y)

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
nonzero_cnt <- colSums(ifelse(out$beta_mat > 0, 1, 0))
plot(out$lambda_seq, nonzero_cnt)

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
res <- microbenchmark::microbenchmark(fitLASSO(X, Y), times = 10)

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
median_time <- median(res$time) / 1e9
# Median Time: 1.29658 sec

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
out <- cvLASSO(X, Y, n_lambda = 30)

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
plot(
  out$lambda_seq,
  out$cvm,
  type = "b",
  pch = 20,
  col = "red"
)
arrows(
  out$lambda_seq,
  out$cvm + out$cvse,
  out$lambda_seq,
  out$cvm - out$cvse,
  angle = 90,
  code = 3,
  length = 0.02,
  col = "gray"
)
