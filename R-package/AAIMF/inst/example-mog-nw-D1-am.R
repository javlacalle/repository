
# NOTE in the amortization approach (linear or NN)
# the estimated round(rprobs, 6) are either 0, 1, see why.

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

# NOTE results below printed for this version
#am_feats <- cbind(1, dgp$X) #dgp$X^2
#am_lst <- list(features = am_feats, fn = NULL, args_fn = list(), finit = NULL, init_args = list())

am_feats <- cbind(1, dgp$X) #dgp$X^2
am_lst <- list(features = am_feats, fn = 'forward_network', args_fn = list(),
               finit = 'init_network', init_args = list(input_dim = NCOL(am_feats)))

model <- make_model(family = 'mog_normal_wishart',
  dgp$X, dgp$K,
  priors = list(),
  amortization = am_lst)

# Initialize the parameters of the model.

model1 <- mog_nw_init(model)

model2 <- mog_init(model,
  seed = 123321,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  #args = list(),
  debug = TRUE)

model <- model2

# ELBO

mog_nw_elbo(model, debug = TRUE)$elbo
# [1] -1059.249

elbo <- mog_elbo(model, do_rprobs = FALSE, debug = TRUE)
elbo$elbo
# [1] -1059.249

# NOTE do_rprobs has no effect,
# rprobs already consistent with etas and pars consistent with rprobs
mog_elbo(model, do_rprobs = TRUE, debug =FALSE)
# [1] -1059.249

# CAVI

# model structure; required to vectorize pars for optim() and unvectorize to build the model
# used internally by mog_cavi_amortized(), actually operates only on model$etas
# model_struct <- get_pars_structure(model$pars)

#model$pars <- mog_nw_upd(model)

#model$rprobs <- mog_nw_upd_resp(model, debug = TRUE)

fit <- mog_cavi_amortized(model,
  optim_args = list(method = "BFGS", hessian = FALSE),
  debug = FALSE)

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
Converged: 0
Iterations: 76 70

Posteriors summary (means)
----------------------------------------------------------------

Mixture probabilities:
[1] 0.3494 0.6506

Component means, m_k:
        [,1]
[1,] -1.9335
[2,]  2.0308

Component plug-in covariances solve(E[Lambda_k])
       [,1]
[1,] 0.4685
       [,1]
[1,] 0.7717

Component mean covariances E[Lambda_k^{-1}]
       [,1]
[1,] 0.4739
       [,1]
[1,] 0.7765

Posterior mean covariance of component means E[Lambda_k^{-1}] / beta_k
       [,1]
[1,] 0.0027
       [,1]
[1,] 0.0024

Final ELBO: -937.0528

Amortization function parameters.
----------------------------------------------------------------
        [,1]
[1,]  1.6487
[2,] 43.1355
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
   2   3 326

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
 Min.   :-1.8881146
 1st Qu.:-0.5062837
 Median :-0.0031333
 Mean   : 0.0001946
 3rd Qu.: 0.5272789
 Max.   : 2.3767407
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
 Min.   :-2.042414
 1st Qu.:-0.506284
 Median :-0.003133
 Mean   :-0.003139
 3rd Qu.: 0.527279
 Max.   : 2.376741
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

# Maux <- matrix(0.0, nrow = dgp$N, ncol = dgp$K)
# for (k in seq_len(dgp$K-1))
#   Maux[,k+1] <- forward_network(X = am_feats, etas = fitted_pars$etas[[k]], return_cache = FALSE)$Y
#
# rprobs <- softmax_rows(Maux)
#rowSums(rprobs)

# the same as model$amortization
fit_am <- fit$model$amortization

if (is.matrix(model$etas))
{

  if (is.null(fit_am$fn))
    fn <- function(X, etas) X %*% etas

  x_grid <- matrix(seq(min(dgp$X), max(dgp$X), length.out = 200), ncol = 1)
  x_grid <- cbind(1, x_grid)
  scores <- do.call(fn, c(list(X = x_grid, etas = cbind(0, fit$model$etas)), fit_am$args_fn))

  plot(x_grid[,2], scores[,2], type = "l")

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

#

# #png("nn-shape.png")
# op <- par(mfrow=c(3,1), mar=c(2,2,2,2))
# plot(x, rprobs[,1])
# plot(x, rprobs[,2])
# plot(x, rprobs[,3])
# par(op)
# #dev.off()

x_grid <- matrix(seq(min(x), max(x), length.out = 200), ncol = 1)
scores <- matrix(0, nrow = nrow(x_grid), ncol = dgp$K)
for (k in seq_len(dgp$K - 1))
  scores[,k+1] <- forward_network(X = x_grid, etas = fitted_pars$etas[[k]], return_cache = FALSE)$Y

#png("nn-shape.png")
op <- par(mfrow=c(2,1), mar=c(2,2,2,2))
plot(x_grid[,1], scores[,2], type = "l")
abline(lm(scores[,2] ~ x_grid[,1]), col = 2)

plot(x_grid[,1], scores[,3], type = "l")
abline(lm(scores[,3] ~ x_grid[,1]), col = 2)
par(op)
#dev.off()

