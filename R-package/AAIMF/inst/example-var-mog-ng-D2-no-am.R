## ============================================================
## Example: CAVI fit for a Generalized VAR model with
## mixture-of-Gaussian innovations and Normal-Gamma priors
## ============================================================

# FIXME
# Check mog_init()
# if (!is.null(ipars$pars$rprobs))
#    {
#      # FIXME check this for AR and no-AR models. Would this block always apply?
#      rprobs <- ipars$pars$rprobs
# }

library(AAIMF)

set.seed(123)

## ------------------------------------------------------------
## 1. Simulate data from a Generalized VAR model
## ------------------------------------------------------------

true_model <- list(
  A = list(
    matrix(c(
        0.55, 0.10,
        0.05, 0.35), nrow = 2, byrow = TRUE),
    matrix(c(
        -0.25, 0.00,
         0.00, -0.10), nrow = 2, byrow = TRUE)),

  means = rbind(
    c(0.0, 0.0),
    c(2.5, 2.0)),

  Sigmas = list(
    matrix(c(
        0.30, 0.08,
        0.08, 0.25), nrow = 2, byrow = TRUE),
    matrix(c(
        0.90, 0.25,
        0.25, 0.70), nrow = 2, byrow = TRUE)),

  probs = c(0.70, 0.30)
)

dgp <- var_mog_sim(n = 800, model = true_model, n.start = 300, seed = 123)

# X <- dgp$X
# Xlags <- dgp$Xlags
# K <- dgp$K
# D <- dgp$D
# p <- dgp$p

cat("\nData dimensions\n")
cat("raw_X:", dim(dgp$raw_X), "\n")
cat("X:", dim(dgp$X), "\n")
cat("Xlags:", dim(dgp$Xlags), "\n")
cat("K:", dgp$K, " D:", dgp$D, " p:", dgp$p, "\n")

cat("\nTrue state counts\n")
print(table(dgp$states))

"
Data dimensions
raw_X: 800 2
X: 798 2
Xlags: 798 4
K: 2  D: 2  p: 2

True state counts

  1   2
561 237
"

# build model

#am_feats <- cbind(1, dgp$X) #dgp$X^2
#am_lst <- list(features = am_feats, fn = NULL, args_fn = list(), finit = NULL, init_args = list())

model <- make_model(family = 'var_mog_normal_gamma',
  dgp$X, dgp$K, p = dgp$p, Xlags = dgp$Xlags,
  priors = list())
  #, amortization = am_lst

# initialization

model1 <- var_mog_ng_init(model,
  init_method = "residual_kmeans",
  kmeans_nstart = 10, kmeans_prob = 0.95,
  ridge_B = 1e-6, ridge_var = 1e-6,
  seed = 123321,
  update_resp = FALSE)

model2 <- mog_init(model,
  #seed = 123321,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  args = list(seed = 125, ridge_B = 1e-6, ridge_var = 1e-8),
  debug = TRUE)

model <- model2

# ELBO

var_mog_ng_elbo(model, TRUE)$elbo
# [1] -2209.226

elbo <- var_mog_ng_elbo(model, TRUE)
elbo$elbo
# [1] -2209.226

elbo2 <- mog_elbo(model, do_rprobs = TRUE, debug = TRUE)
elbo2$elbo
# [1] -2209.226

# CAVI

#model$pars <- var_mog_ng_upd(model)

#model$rprobs <- var_mog_ng_upd_resp(model, debug = TRUE)

#fit <- mog_cavi_amortized(model,
#  optim_args = list(method = "BFGS", hessian = FALSE),
#  debug = FALSE)
#fit$posteriors

fit <- mog_cavi(model,
  maxiter = 250, tol = 1e-8,
  verbose = TRUE,
  debug = FALSE)

"
iter =    1 | ELBO = -1978.051419
iter =    2 | ELBO = -1953.597987
iter =    3 | ELBO = -1949.049032
iter =    4 | ELBO = -1947.321694
iter =    5 | ELBO = -1946.648329
iter =    6 | ELBO = -1946.388828
iter =    7 | ELBO = -1946.289390
iter =    8 | ELBO = -1946.251460
iter =    9 | ELBO = -1946.237072
iter =   10 | ELBO = -1946.231653
iter =   11 | ELBO = -1946.229628
iter =   12 | ELBO = -1946.228879
iter =   13 | ELBO = -1946.228604
iter =   14 | ELBO = -1946.228504
iter =   15 | ELBO = -1946.228468
iter =   16 | ELBO = -1946.228455
"

fit

"
Variational Bayes fit for VAR mixture-of-Gaussians-noise model
==============================================================

Model dimensions
----------------------------------------------------------------
Number of observations used: 798
AR order p: 2
Number of mixture components K: 2

Convergence
----------------------------------------------------------------
Converged: TRUE
Iterations: 16
Initial ELBO: -2209.226
Final ELBO: -1946.228

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities: 0.2889 0.7111

Innovation means, m_k:
         x1      x2
[1,] 2.5466  2.0218
[2,] 0.0043 -0.0098

Component expected precisions E[lambda_kd]:
        x1     x2
[1,] 1.164 1.4362
[2,] 3.465 4.3237

Innovation plug-in variances 1 / E[lambda_kd]
[1] 0.8591
[1] 0.2886
[1] 0.6963
[1] 0.2313

Innovation mean variances E[1 / lambda_kd]
[1] 0.8665
[1] 0.2896
[1] 0.7023
[1] 0.2321

VAR coefficient matrices (Bbar converted to A matrices)
       [,1]   [,2]
[1,] 0.5304 0.1275
[2,] 0.0409 0.3523
        [,1]    [,2]
[1,] -0.2266 -0.0074
[2,] -0.0036 -0.0621
"

# True values

cat("\nMixture probabilities\n")
print(true_model$probs)

cat("\nComponent innovation means\n")
print(true_model$means)

cat("\nInnovation covariance matrices\n")
print(true_model$Sigmas)

cat("\nVAR coefficient matrices\n")
print(true_model$A)

"
Mixture probabilities
[1] 0.7 0.3

Component innovation means
     [,1] [,2]
[1,]  0.0    0
[2,]  2.5    2

Innovation covariance matrices
[[1]]
     [,1] [,2]
[1,] 0.30 0.08
[2,] 0.08 0.25

[[2]]
     [,1] [,2]
[1,] 0.90 0.25
[2,] 0.25 0.70


VAR coefficient matrices
[[1]]
     [,1] [,2]
[1,] 0.55 0.10
[2,] 0.05 0.35

[[2]]
      [,1] [,2]
[1,] -0.25  0.0
[2,]  0.00 -0.1
"

# state classification diagnostics

z_hat <- max.col(fit$model$rprobs)

cat("\nClassification table, without label correction\n")
print(table(true = dgp$states, estimated = z_hat))

if (dgp$K == 2L)
{
  acc_identity <- mean(z_hat == dgp$states)
  acc_swapped <- mean((3L - z_hat) == dgp$states)

  cat("\nBest label-adjusted accuracy:", max(acc_identity, acc_swapped), "\n")
}

"
    estimated
true   1   2
   1   6 555
   2 222  15

Best label-adjusted accuracy: 0.9736842
"

# ELBO plot

plot(fit$elbo_path,
  type = "l", xlab = "Iteration", ylab = "ELBO",
  main = "Generalized VAR MoG Normal-Gamma CAVI fit")

# residual summary under posterior mean Bbar
# FIXME (see include means to get zero-mean residuals)
# see example-mog-nw-D1-no-am 'residual summary under hard component assignment'

Xhat <- dgp$Xlags %*% fit$model$pars$Bbar
Ehat <- dgp$X - Xhat

cat("\nResidual summary under posterior mean Bbar\n")
print(summary(Ehat))

"
Residual summary under posterior mean Bbar
       x1                x2
 Min.   :-1.5162   Min.   :-1.8812
 1st Qu.:-0.1959   1st Qu.:-0.2028
 Median : 0.3027   Median : 0.2376
 Mean   : 0.7375   Mean   : 0.5761
 3rd Qu.: 1.5405   3rd Qu.: 1.1402
 Max.   : 5.3243   Max.   : 4.1541
"

# inspect posterior responsibilities over time

matplot(fit$model$rprobs, type = "l", lty = 1,
  xlab = "Time index after lag trimming", ylab = "Responsibility",
  main = "Posterior responsibilities")

legend( "topright", legend = paste0("component ", seq_len(model$K)),
  col = seq_len(model$K), lty = 1, bty = "n")

# zoom in into the first observations
plot(ts(fit$model$rprobs[1:100,]),type="s")
