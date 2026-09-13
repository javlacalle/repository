library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(-2, 2), Sigmas = list(matrix(0.4, 1, 1), matrix(0.8, 1, 1)),
  seed = 123)
#X <- dgp$X
#print(dgp)

cat("\nData dimensions\n")
cat("raw_X:", dim(dgp$raw_X), "\n")
cat("X:", dim(dgp$X), "\n")
cat("K:", dgp$K, " D:", dgp$D, "\n")

cat("\nTrue state counts\n")
print(table(dgp$states))

"
Data dimensions
raw_X:
X: 500 1
K: 2  D: 1

True state counts

  1   2
171 329
"

# build model

model <- make_model(family = 'mog_normal_wishart',
  dgp$X, dgp$K,
  priors = list())

# initialization

model1 <- mog_nw_init(model)

model2 <- mog_init(model,
  seed = 123321,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  #args = list(),
  debug = TRUE)

model <- model2

# ELBO

mog_nw_elbo(model, TRUE)$elbo
# [1] -1031.255

elbo <- mog_elbo(model,
  amortization = list(features = NULL, fn = NULL, args = list()),
  debug = TRUE)
elbo$elbo
# [1] -1031.255

# CAVI

#model$pars <- mog_nw_upd(model)

#model$rprobs <- mog_nw_upd_resp(model, debug = TRUE)

fit <- mog_cavi(model,
  maxiter = 250, tol = 1e-8,
  verbose = TRUE,
  debug = FALSE)

"
iter =    1 | ELBO = -964.585131
iter =    2 | ELBO = -944.086308
iter =    3 | ELBO = -937.367517
iter =    4 | ELBO = -935.006528
iter =    5 | ELBO = -934.450343
iter =    6 | ELBO = -934.346453
iter =    7 | ELBO = -934.328261
iter =    8 | ELBO = -934.325122
iter =    9 | ELBO = -934.324583
iter =   10 | ELBO = -934.324490
iter =   11 | ELBO = -934.324475
iter =   12 | ELBO = -934.324472
"

fit

"
Variational Bayes fit for mixture-of-Gaussians model
====================================================

Model dimensions
----------------------------------------------------------------
Number of observations: 500
Number of variables D: 1
Number of mixture components K: 2

Convergence
----------------------------------------------------------------
Converged: TRUE
Iterations: 12
Initial ELBO: -1031.255
Final ELBO: -934.3245

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities:
[1] 0.6548 0.3452

Component means, m_k:
        [,1]
[1,]  2.0140
[2,] -1.9502

Component plug-in covariances solve(E[Lambda_k])
       [,1]
[1,] 0.8103
       [,1]
[1,] 0.4515

Component mean covariances E[Lambda_k^{-1}]
       [,1]
[1,] 0.8153
       [,1]
[1,] 0.4567

Posterior mean covariance of component means E[Lambda_k^{-1}] / beta_k
       [,1]
[1,] 0.0025
       [,1]
[1,] 0.0026
"

# True values

cat("\nMixture probabilities\n")
print(dgp$pis)

cat("\nComponent innovation means\n")
print(dgp$mus)

cat("\nInnovation covariance matrices\n")
print(dgp$Sigmas)

"
Mixture probabilities
[1] 0.35 0.65

Component innovation means
     [,1]
[1,]   -2
[2,]    2

Innovation covariance matrices
[[1]]
     [,1]
[1,]  0.4

[[2]]
     [,1]
[1,]  0.8
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
   1   0 171
   2 326   3

Best label-adjusted accuracy: 0.994
"

# ELBO plot
# NOTE not available for neural-network amortization function

if (!is.null(fit$elbo_path))
  plot(fit$elbo_path,
    type = "l", xlab = "Iteration", ylab = "ELBO",
    main = "MoG Normal-Wishart CAVI fit")

# residual summary under posterior mean Bbar

mu_hat_n <- fit$model$rprobs %*% fit$model$pars$mu_ms
Ehat <- dgp$X - mu_hat_n

cat("\nResidual summary under posterior responsibility-weighted component means\n")
print(summary(Ehat))

"
Residual summary under posterior responsibility-weighted component means
       V1
 Min.   :-1.7215060
 1st Qu.:-0.4894902
 Median : 0.0082175
 Mean   : 0.0001276
 3rd Qu.: 0.5156499
 Max.   : 2.3935879
"

# residual summary under hard component assignment

z_hat <- max.col(fit$model$rprobs)
mu_hat_hard <- fit$model$pars$mu_ms[z_hat,,drop=FALSE]
Ehat_hard <- dgp$X - mu_hat_hard

cat("\nResidual summary under hard posterior component assignment\n")
print(summary(Ehat_hard))

"
Residual summary under hard posterior component assignment
       V1
 Min.   :-2.02557
 1st Qu.:-0.48949
 Median : 0.01362
 Mean   : 0.01364
 3rd Qu.: 0.54413
 Max.   : 2.39359
"

# inspect posterior responsibilities over time

matplot(fit$model$rprobs, type = "l", lty = 1,
  xlab = "Time index after lag trimming", ylab = "Responsibility",
  main = "Posterior responsibilities")

legend( "topright", legend = paste0("component ", seq_len(model$K)),
  col = seq_len(model$K), lty = 1, bty = "n")

# zoom-in on the first observations

plot(ts(fit$model$rprobs[1:100,]),type="s")
