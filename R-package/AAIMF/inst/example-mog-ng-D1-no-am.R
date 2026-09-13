library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(-2, 2), Sigmas = list(matrix(0.4, 1, 1), matrix(0.8, 1, 1)),
  seed = 123)

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

model <- make_model(family = 'mog_normal_gamma',
  dgp$X, dgp$K,
  priors = list())

# initialization

model1 <- mog_nw_init(model)

model2 <- mog_init(model,
  seed = 123321,
  amortization = list(features = NULL, fn = NULL, args_fn = list(), finit = NULL, init_args = list()),
  kmeans_nstart = 10, kmeans_prob = 0.95,
  #args = list(),
  debug = TRUE)

model <- model2

# ELBO

mog_ng_elbo(model, TRUE)$elbo
[1] -1030.496

elbo <- mog_elbo(model,
  amortization = list(features = NULL, fn = NULL, args = list()),
  debug = TRUE)
elbo$elbo
# [1] -1030.496

# CAVI

#model$pars <- mog_ng_upd(model)

#model$rprobs <- mog_ng_upd_resp(model, debug = TRUE)

fit <- mog_cavi(model,
  maxiter = 250, tol = 1e-8,
  verbose = TRUE,
  debug = FALSE)

"
iter =    1 | ELBO = -963.667791
iter =    2 | ELBO = -943.267372
iter =    3 | ELBO = -936.659651
iter =    4 | ELBO = -934.366565
iter =    5 | ELBO = -933.830567
iter =    6 | ELBO = -933.730995
iter =    7 | ELBO = -933.713669
iter =    8 | ELBO = -933.710701
iter =    9 | ELBO = -933.710195
iter =   10 | ELBO = -933.710108
iter =   11 | ELBO = -933.710094
iter =   12 | ELBO = -933.710091
"

fit

"
Variational Bayes fit for mixture-of-Gaussians model
Normal-Gamma prior with diagonal component precisions
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
Initial ELBO: -1030.496
Final ELBO: -933.7101

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities:
[1] 0.6547 0.3453

Component means, m_k:
        [,1]
[1,]  2.0143
[2,] -1.9496

Component expected precisions E[lambda_kd]:
       [,1]
[1,] 1.2343
[2,] 2.1954

Component plug-in variances 1 / E[lambda_kd]:
       [,1]
[1,] 0.8102
[2,] 0.4555

Component mean variances E[1 / lambda_kd]:
       [,1]
[1,] 0.8151
[2,] 0.4608

Posterior mean variances of component means E[1 / lambda_kd] / beta_k:
       [,1]
[1,] 0.0025
[2,] 0.0027
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

plot(fit$elbo_path,
  type = "l", xlab = "Iteration", ylab = "ELBO",
  main = "MoG Normal-Gamma CAVI fit")

# residual summary under posterior mean Bbar

mu_hat_n <- fit$model$rprobs %*% fit$model$pars$mu_ms
Ehat <- dgp$X - mu_hat_n

cat("\nResidual summary under posterior responsibility-weighted component means\n")
print(summary(Ehat))

"
Residual summary under posterior responsibility-weighted component means
       V1
 Min.   :-1.7220639
 1st Qu.:-0.4899073
 Median : 0.0105406
 Mean   : 0.0001294
 3rd Qu.: 0.5152350
 Max.   : 2.3932171
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
 Min.   :-2.02594
 1st Qu.:-0.48990
 Median : 0.01315
 Mean   : 0.01320
 3rd Qu.: 0.54376
 Max.   : 2.39322
"

# inspect posterior responsibilities over time

matplot(fit$model$rprobs, type = "l", lty = 1,
  xlab = "Time index after lag trimming", ylab = "Responsibility",
  main = "Posterior responsibilities")

legend( "topright", legend = paste0("component ", seq_len(model$K)),
  col = seq_len(model$K), lty = 1, bty = "n")

#plot(ts(fit$model$rprobs[1:100,]),type="s")
