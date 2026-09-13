
# NOTE in the amortization approach (linear or NN)
# the estimated round(rprobs, 6) are either 0, 1, see why.

library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(c(-2, -1), c(2, 1)),
  Sigmas = list(diag(c(0.4, 0.8)), matrix(c(0.8, 0.3, 0.3, 0.6), 2, 2)),
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
X: 500 2
K: 2  D: 2

True state counts

  1   2
171 329
"

# build model

am_feats <- cbind(1, dgp$X) #dgp$X^2
am_lst <- list(features = am_feats, fn = NULL, args_fn = list(), finit = NULL, init_args = list())

#am_feats <- cbind(1, dgp$X) #dgp$X^2
#am_lst <- list(features = am_feats, fn = 'forward_network', args_fn = list(),
#               finit = 'init_network', init_args = list(input_dim = NCOL(am_feats)))

model <- make_model(family = 'mog_normal_gamma',
  dgp$X, dgp$K,
  priors = list(),
  amortization = am_lst)

# Initialize the parameters of the model.

model1 <- mog_ng_init(model)

model2 <- mog_init(model,
  seed = 123321,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  #args = list(),
  debug = TRUE)

model <- model2

# ELBO

mog_ng_elbo(model, debug = TRUE)$elbo
# [1] -1730.844

elbo <- mog_elbo(model, do_rprobs = FALSE, debug = TRUE)
elbo$elbo
# [1] -1730.844

# NOTE do_rprobs has no effect,
# rprobs already consistent with etas and pars consistent with rprobs
mog_elbo(model, do_rprobs = TRUE, debug =FALSE)
# [1] -1730.844

# CAVI

# model structure; required to vectorize pars for optim() and unvectorize to build the model
# used internally by mog_cavi_amortized(), actually operates only on model$etas
# model_struct <- get_pars_structure(model$pars)

#model$pars <- mog_ng_upd(model)

#model$rprobs <- mog_ng_upd_resp(model, debug = TRUE)

fit <- mog_cavi_amortized(model,
  optim_args = list(method = "BFGS", hessian = FALSE),
  debug = FALSE)

fit

"
Variational Bayes fit for mixture-of-Gaussians model
Normal-Gamma prior with diagonal component precisions
====================================================

Model dimensions
----------------------------------------------------------------
Number of observations: 500
Number of variables D: 2
Number of mixture components K: 2

Convergence
----------------------------------------------------------------
Converged: 0
Iterations: 4 2

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities:
[1] 0.3526 0.6474

Component means, m_k:
        [,1]    [,2]
[1,] -1.9088 -0.9400
[2,]  2.0527  1.0252

Component expected precisions E[lambda_kd]:
       [,1]   [,2]
[1,] 1.8941 1.3612
[2,] 1.1983 1.8047

Component plug-in variances 1 / E[lambda_kd]:
       [,1]   [,2]
[1,] 0.5279 0.7346
[2,] 0.8345 0.5541

Component mean variances E[1 / lambda_kd]:
       [,1]   [,2]
[1,] 0.5339 0.7430
[2,] 0.8397 0.5575

Posterior mean variances of component means E[1 / lambda_kd] / beta_k:
       [,1]   [,2]
[1,] 0.0030 0.0042
[2,] 0.0026 0.0017

Final ELBO: -1557.056

Amortization function parameters.
----------------------------------------------------------------
         [,1]
[1,] -26.4578
[2,] 176.2959
[3,] 102.0389
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
   1 171   0
   2   5 324

Best label-adjusted accuracy: 0.99
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
       V1                   V2
 Min.   :-2.0299118   Min.   :-2.5731407
 1st Qu.:-0.5290607   1st Qu.:-0.5338996
 Median :-0.0114342   Median :-0.0243849
 Mean   : 0.0002877   Mean   : 0.0001704
 3rd Qu.: 0.5155864   3rd Qu.: 0.5266319
 Max.   : 2.3548807   Max.   : 2.3357736

"

# residual summary under hard component assignment

z_hat <- max.col(fit$model$rprobs)
mu_hat_hard <- fit$model$pars$mu_ms[z_hat,,drop=FALSE]
Ehat_hard <- dgp$X - mu_hat_hard

cat("\nResidual summary under hard posterior component assignment\n")
print(summary(Ehat_hard))

"
       V1                   V2
 Min.   :-2.0299118   Min.   :-2.5731407
 1st Qu.:-0.5290607   1st Qu.:-0.5338996
 Median :-0.0114342   Median :-0.0243849
 Mean   : 0.0002876   Mean   : 0.0001704
 3rd Qu.: 0.5155864   3rd Qu.: 0.5266319
 Max.   : 2.3548807   Max.   : 2.3357736
"

# inspect posterior responsibilities over time

matplot(fit$model$rprobs, type = "l", lty = 1,
  xlab = "Time index after lag trimming", ylab = "Responsibility",
  main = "Posterior responsibilities")

legend( "topright", legend = paste0("component ", seq_len(model$K)),
  col = seq_len(model$K), lty = 1, bty = "n")

# zoom-in on the first observations

plot(ts(fit$model$rprobs[1:100,]),type="s")

# Inspect the function learned by the neural-network.
# TODO in the example with NN amortization instead of linear function

# the same as model$amortization
fit_am <- fit$model$amortization

if (is.matrix(model$etas))
{
  if (is.null(fit_am$fn))
    fn <- function(X, etas) X %*% etas

  x1_grid <- matrix(seq(min(dgp$X[,1]), max(dgp$X[,1]), length.out = 200), ncol = 1)
  x2_grid <- matrix(seq(min(dgp$X[,2]), max(dgp$X[,2]), length.out = 200), ncol = 1)
  x_grid <- cbind(1, x1_grid, x2_grid)
  scores <- do.call(fn, c(list(X = x_grid, etas = cbind(0, fit$model$etas)), fit_am$args_fn))

  plot(x_grid[,2], scores[,2], type = "l")
  plot(x_grid[,3], scores[,2], type = "l")

} else
if (is.list(model$etas))
{
  R <- length(model$etas) + 1L
  scores <- matrix(0.0, nrow(fit_am$features), R)

  x_grid <- matrix(seq(min(dgp$X), max(dgp$X), length.out = 200), ncol = 1)
  x_grid <- cbind(1, x_grid)
  scores <- matrix(0, nrow = nrow(x_grid), ncol = dgp$K)

  for (r in seq_len(R - 1))
    scores[,r+1] <- do.call(fit_am$fn,
      c(list(X = x_grid, etas = model$etas[[r]]), fit_am$args_fn))$Y

  #plot(x_grid[,1], scores[,1], type = "l")
  #plot(x_grid[,1], scores[,2], type = "l")
  plot(x_grid[,2], scores[,1], type = "l")
  plot(x_grid[,2], scores[,2], type = "l")
}
