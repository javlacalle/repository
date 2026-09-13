## ============================================================
## Example: CAVI fit for a Generalized VAR model with
## mixture-of-Gaussian innovations and Normal-Wishart priors
## ============================================================

# FIXME
# Note that the example assumes the updated <code>mog_cavi()</code> accepts <code>init_args</code> and passes them to the family initializer. If your patched generic function uses a slightly different argument name, only the fitting block needs to be adjusted.

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

# build model


#am_feats <- cbind(1, dgp$X) #dgp$X^2
#am_lst <- list(features = am_feats, fn = NULL, args_fn = list(), finit = NULL, init_args = list())


model <- make_model(family = 'var_mog_normal_wishart',
  dgp$X, dgp$K, p = dgp$p, Xlags = dgp$Xlags,
  priors = list())
  #, amortization = am_lst

# initialization

model1 <- var_mog_nw_init(model,
  init_method = "residual_kmeans",
  kmeans_nstart = 10, kmeans_prob = 0.95,
  ridge_B = 1e-6, ridge_cov = 1e-6,
  seed = 123321,
  update_resp = FALSE)

model2 <- mog_init(model,
  seed = 123321,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  args = list(seed = 125),
  debug = TRUE)

model <- model2

# FIXME inspect; Note that the second time the argument model has been modified by var_mog_nw_init()
for (label in names(model))
{
  if (!identical(model1[[label]], model2[[label]]))
    print(label)
}
# [1] "rprobs"
# [1] "pars"

# ELBO

var_mog_nw_elbo(model1, TRUE)$elbo
# [1] -2140.662

elbo <- var_mog_nw_elbo(model2, TRUE)
elbo$elbo
# [1] -2081.806

elbo1 <- mog_elbo(model1,
  amortization = list(features = NULL, fn = NULL, args = list()),
  debug = TRUE)
elbo1$elbo
# [1] -2140.662

elbo2 <- mog_elbo(model, do_rprobs = TRUE, debug = TRUE)
elbo2$elbo
# [1] -2081.806

# CAVI

#model$pars <- var_mog_nw_upd(model)

#model$rprobs <- var_mog_nw_upd_resp(model, debug = TRUE)

#fit <- mog_cavi_amortized(model,
#  optim_args = list(method = "BFGS", hessian = FALSE),
#  debug = FALSE)
#fit$posteriors

fit <- mog_cavi(model,
  maxiter = 250, tol = 1e-8,
  verbose = TRUE,
  debug = FALSE)

"
iter =    1 | ELBO = -1953.175976
iter =    2 | ELBO = -1924.927439
iter =    3 | ELBO = -1919.169635
iter =    4 | ELBO = -1917.835314
iter =    5 | ELBO = -1917.411506
iter =    6 | ELBO = -1917.251036
iter =    7 | ELBO = -1917.186434
iter =    8 | ELBO = -1917.159863
iter =    9 | ELBO = -1917.148844
iter =   10 | ELBO = -1917.144258
iter =   11 | ELBO = -1917.142346
iter =   12 | ELBO = -1917.141549
iter =   13 | ELBO = -1917.141217
iter =   14 | ELBO = -1917.141079
iter =   15 | ELBO = -1917.141021
iter =   16 | ELBO = -1917.140997
iter =   17 | ELBO = -1917.140987
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
Iterations: 17
Initial ELBO: -2004.465
Final ELBO: -1917.141

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities: 0.703 0.297

Innovation means, m_k:
        [,1]    [,2]
[1,] -0.0005 -0.0126
[2,]  2.4958  1.9832

Innovation plug-in covariances solve(E[Lambda_k])
       [,1]   [,2]
[1,] 0.2806 0.0576
[2,] 0.0576 0.2250
       [,1]   [,2]
[1,] 0.9300 0.3182
[2,] 0.3182 0.7319

Innovation mean covariances E[Lambda_k^{-1}]
       [,1]   [,2]
[1,] 0.2821 0.0579
[2,] 0.0579 0.2262
       [,1]   [,2]
[1,] 0.9418 0.3222
[2,] 0.3222 0.7411

VAR coefficient matrices (Bbar converted to A matrices)
       [,1]   [,2]
[1,] 0.5274 0.1306
[2,] 0.0368 0.3559
        [,1]    [,2]
[1,] -0.2225 -0.0144
[2,]  0.0013 -0.0703
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
Classification table, without label correction
    estimated
true   1   2
   1 555   6
   2  14 223

Best label-adjusted accuracy: 0.9749373
"

# ELBO plot

plot(fit$elbo_path,
  type = "l", xlab = "Iteration", ylab = "ELBO",
  main = "Generalized VAR MoG Normal-Wishart CAVI fit")

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
 Min.   :-1.5118   Min.   :-1.8842
 1st Qu.:-0.1983   1st Qu.:-0.1992
 Median : 0.3058   Median : 0.2362
 Mean   : 0.7396   Mean   : 0.5791
 3rd Qu.: 1.5368   3rd Qu.: 1.1480
 Max.   : 5.3243   Max.   : 4.1642
"

# inspect posterior responsibilities over time

matplot(fit$model$rprobs, type = "l", lty = 1,
  xlab = "Time index after lag trimming", ylab = "Responsibility",
  main = "Posterior responsibilities")

legend( "topright", legend = paste0("component ", seq_len(model$K)),
  col = seq_len(model$K), lty = 1, bty = "n")

#plot(ts(fit$model$rprobs[1:100,]),type="s")
