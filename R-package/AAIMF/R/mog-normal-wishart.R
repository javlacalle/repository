## Variational Bayes for the
## Normal-Wishart mixture-of-Gaussians-noise model.


## -------------------------------------------------------------------
## Initialization of parameters.
## -------------------------------------------------------------------

mog_nw_init <- function(model)
{
  # Build the model with make_model() for preliminary checks.

  # In this model, only priors need to be specifically initialized here.
  # The remaining parameters are initialized by the general-purpose function mog_init().

  priors <- model$priors
  K <- model$K
  D <- model$D

  if (length(priors$alpha0) == 1) {
    priors$alpha0 <- rep(priors$alpha0, K)
  } else
  if (length(priors$alpha0) != K)
    stop("'priors$alpha0' must have length 1 or K.")

  if (is.null(priors$m0)) {
    priors$m0 <- rep(0, D)
  } else
  if (length(priors$m0) != D)
    stop("'priors$m0' must have length D.")

  if (is.null(priors$W0)) {
    priors$W0 <- diag(D)
  } else
  if (!all(dim(priors$W0) == c(D, D)))
    stop("'priors$W0' must be a D x D matrix.")

  priors$W0_inv <- solve(priors$W0)

  if (is.null(priors$nu0)) {
    priors$nu0 <- D
  } else
  if (priors$nu0 <= D - 1)
    stop("'priors$nu0' must be greater than D - 1 for a valid Wishart density.")

  # Not needed. This is done by make_model() which ensure proper arguments
  # priors <- modifyList_priorsPars(formals(mog_nw_init)[[1]], priors)

  #list(priors = priors, rprobs = NULL, pars = NULL)
  model$priors <- priors
  model
}


## -------------------------------------------------------------------
## CAVI updates.
## -------------------------------------------------------------------

mog_nw_upd_aux <- function(X, rprobs, colSums_rprobs, K, D)
{
  # Internal function. Assumes callers pass valid arguments.

  eps <- 1e-12

  xbar <- matrix(0, nrow = K, ncol = D)
  S <- vector("list", K)
  S[] <- list(matrix(0, D, D)) # initialize all matrices with zeros

  for (k in seq_len(K))
  {
    if (colSums_rprobs[k] > eps)
    {
      xbar[k,] <- colSums(rprobs[,k] * X) / colSums_rprobs[k]
      #Xcentered <- sweep(X, MARGIN = 2, STATS = xbar[k,], FUN = "-", check.margin = FALSE)
      Xcentered <- sweep(X, 2, xbar[k,], "-", FALSE)
      S[[k]] <- crossprod(Xcentered, rprobs[,k] * Xcentered) / colSums_rprobs[k]
      #crossprod(Xcentered, sweep(Xcentered, 1L, rprobs[, k], "*")) / colSums_rprobs[k]

    } # else (initialized with zeros)
  }

  list(xbar = xbar, S = S)
}

mog_nw_upd <- function(model)
{
  # Internal function. Assumes callers pass valid arguments.

  X <- model$X
  priors <- model$priors
  rprobs <- model$rprobs
  K <- model$K
  D <- model$D

  # FIXME See pass 'Lambda_W', ... as argument and modify, instead of creating here all the containers every time.

  #eps <- 1e-12
  colSums_rprobs <- colSums(rprobs)
  aux <- mog_nw_upd_aux(X, rprobs, colSums_rprobs, K, D)

  pi_alphas <- priors$alpha0 + colSums_rprobs

  mu_betas <- priors$beta0 + colSums_rprobs
  mu_ms <- (priors$beta0 * matrix(priors$m0, K, D, byrow=TRUE) +
    matrix(colSums_rprobs, K, D, byrow=FALSE) * aux$xbar) / mu_betas

  Winv_weigths <- priors$beta0 * colSums_rprobs / mu_betas
  Lambda_nus <- priors$nu0 + colSums_rprobs

  #mu_ms <- matrix(0, nrow = K, ncol = D)
  # Old: The names are required by modifyList() within mog2_cavi().
  Lambda_Ws <- setNames(vector("list", K), paste0("k", seq_len(K)))

  for (k in seq_len(K))
  {
    # Do not consider here the cases
    # if (colSums_rprobs[k] > eps) { ... } else { ... }.
    # mog2_update_auxiliary() sets colSums_rprobs[k] to 0 when < eps;
    # so it may be better to keep here multiplication by zero rather than
    # adding if-else within this loop (also 'm' can be vectorized outside the loop).

    # FIXME debug
    #m[k,] <- (priors$beta0 * priors$m0 + colSums_rprobs[k] * aux$xbar[k,]) / mu_betas[k]

    xbar_centered <- sweep(aux$xbar[k,,drop=FALSE], 2L, priors$m0, "-", FALSE)
    W_inv <- priors$W0_inv + colSums_rprobs[k] * aux$S[[k]] +
      Winv_weigths[k] * crossprod(xbar_centered)
    Lambda_Ws[[k]] <- solve(W_inv)
  }

  list(pi_alphas = pi_alphas, mu_ms = mu_ms, mu_betas = mu_betas,
       Lambda_Ws = Lambda_Ws, Lambda_nus = Lambda_nus)
}

mog_nw_upd_resp <- function(model, debug = TRUE)
{
  # Internal function. Assumes callers pass valid arguments.
  # Recall: Amortization is not addressed by CAVI.

  X <- model$X
  pars <- model$pars
  K <- model$K
  N <- model$N
  D <- model$D

  Dlog2 <- D * log(2)
  Eq_log_pi <- digamma(pars$pi_alphas) - digamma(sum(pars$pi_alphas))
  Eq_log_Lambda <- vapply(seq_len(K), function(k)
    expected_log_det_wishart(pars$Lambda_Ws[[k]], pars$Lambda_nus[[k]], D, Dlog2),
    numeric(1))
  hELmDlog2pi <- 0.5 * (Eq_log_Lambda - D * log(2 * pi))

  scores <- matrix(0, N, K)

  # FIXME Same Qnk as in ELBO, see reuse as separate function.

  for (k in seq_len(K))
  {
    Qnk <- D / pars$mu_betas[k] +
      pars$Lambda_nus[k] * quadratic_byrows(X, pars$mu_ms[k,], pars$Lambda_Ws[[k]])

    scores[,k] <- Eq_log_pi[k] + hELmDlog2pi[k] - 0.5 * Qnk
  }

  rprobs <- softmax_rows(scores)

  if (debug && !isTRUE(all.equal(rowSums(rprobs), rep(1, N))))
    stop("Responsibilities do not sum to one.")

  rprobs
}


## -------------------------------------------------------------------
## Evidence lower bound (ELBO) for the MoG model.
## Analytical ELBO.
## -------------------------------------------------------------------

quadratic_byrows <- function(X, center, A)
{
  # Based on stats::mahalanobis().
  #Xcentered <- sweep(x = X, MARGIN = 2L, STATS = center, FUN = "-", check.margin = FALSE)
  X_centered <- sweep(X, 2L, center, "-", FALSE)
  rowSums((X_centered %*% A) * X_centered)
}

logdet_spd_mat <- function(A)
{
  # Logarithm of determinant of symmetric positive definite matrix.

  # For matrices of order > 50 x 50 the following is faster.
  # (The larger the matrix the larger the difference in timings)
  #R <- tryCatch(chol(A), error = function(e) stop("Matrix is not positive definite: ", e$message))
  #2 * sum(log(diag(R)))

  # For matrices of order < 50 x 50 this seems faster.
  # (For orders close to 50 x 50 similar timings.)

  detA <- determinant(A, logarithm = TRUE)
  if (detA$sign <= 0)
    stop("Matrix is not positive definite.")
  as.numeric(detA$modulus)
}

log_wishart_constant <- function(W, nu, hDnulog2, qDDm1logpi, sum_hlg_np1mseqD)
{
   # Logarithm of the Wishart constant in a
   # Wishart distribution with scale W and degrees of freedom nu.

# print(-0.5 * nu * logdet_spd_mat(W))
# print(hDnulog2)
# print(qDDm1logpi)
# print(sum_hlg_np1mseqD)

  -0.5 * nu * logdet_spd_mat(W) - hDnulog2 - qDDm1logpi - sum_hlg_np1mseqD
}

expected_log_det_wishart <- function(W, nu, D, Dlog2)
{
#print(logdet_spd_mat(W))
  sum(digamma((nu + 1 - seq_len(D)) / 2)) + Dlog2 + logdet_spd_mat(W)
}

entropy_wishart <- function(W, nu, D, hD, Dlog2, hDnulog2, qDDm1logpi, sum_hlg_np1mseqD, hnumDm1)
{
# print(-log_wishart_constant(W, nu, hDnulog2, qDDm1logpi, sum_hlg_np1mseqD))
# print(-hnumDm1)
# print(expected_log_det_wishart(W, nu, D, Dlog2))
# print(hD * nu)

  -log_wishart_constant(W, nu, hDnulog2, qDDm1logpi, sum_hlg_np1mseqD) -
    hnumDm1 * expected_log_det_wishart(W, nu, D, Dlog2) + hD * nu
}

entropy_dirichlet <- function(alphas, K)
{
  a <- sum(alphas)
  sum(lgamma(alphas)) - lgamma(a) +
  (a - K) * digamma(a) - sum((alphas - 1) * digamma(alphas))
}

mog_nw_elbo <- function(model, debug = TRUE)
{
  X <- model$X
  priors <- model$priors
  rprobs <- model$rprobs
  pars <- model$pars
  K <- model$K
  D <- model$D

  # Constants.

  #colSums_rprobs <- colSums(pars$rprobs)

  log2pi <- log(2 * pi)
  hD <- 0.5 * D
  hDlog2pi <- hD * log2pi
  log_beta0 <- hD * log(priors$beta0)
  hbeta0 <- 0.5 * priors$beta0
  hDlogbeta0 <- hD * log(priors$beta0)
  hnu0mdp1 <- 0.5 * (priors$nu0 - D - 1)
  hnus <- 0.5 * pars$Lambda_nus
  hDnu0log2 <- hD * priors$nu0 * log(2)
  hDnulog2 <- hD * pars$Lambda_nus * log(2)
  qDDm1logpi <- 0.25 * D * (D - 1) * log(pi)
  sum_hlg_n0p1mseqD <- sum(lgamma((priors$nu0 + 1 - seq_len(D)) / 2))
  sum_hlg_np1mseqD <- vapply(X = seq_len(K), FUN = function(k)
    sum(lgamma((pars$Lambda_nus[k] + 1 - seq_len(D)) / 2)), numeric(1))
  hD_1plog2pi <- hD * (1 + log2pi)
  hnumDm1 <- 0.5 * (pars$Lambda_nus - D - 1)
  Dlog2 <- D * log(2)

  # Expected log-priors.

  Eq_log_pi <- digamma(pars$pi_alphas) - digamma(sum(pars$pi_alphas))
  Eq_log_Lambda <- vapply(seq_len(K), function(k)
    expected_log_det_wishart(pars$Lambda_Ws[[k]], pars$Lambda_nus[[k]], D, Dlog2),
    numeric(1))

#print_named(Eq_log_pi, "res-v2.txt", FALSE)
#print_named(Eq_log_Lambda, "res-v2.txt")

  # Expected log likelihood E_q[log p(X | Z, mu, Lambda)].

  Eq_log_lik <- 0
  for (k in seq_len(K))
  {
    Qnk <- D / pars$mu_betas[k] +
      pars$Lambda_nus[k] * quadratic_byrows(X, pars$mu_ms[k,], pars$Lambda_Ws[[k]])
    Eq_log_lik <- Eq_log_lik +
      sum(rprobs[,k] * (0.5 * Eq_log_Lambda[k] - hDlog2pi - 0.5 * Qnk))
  }

#print_named(Eq_log_lik, "res-v2.txt")

  # Expected log p(Z | pi).

  Eq_log_Z <- sum(sweep(rprobs, 2L, Eq_log_pi, "*", FALSE))
  #sum(sweep(x = rprobs, MARGIN = 2L, STATS = Eq_log_pi, FUN = "*", check.margin = FALSE))

#print_named(Eq_log_Z, "res-v2.txt")

  # Expected log p(pi).

  Eq_log_prior_pi <- lgamma(sum(priors$alpha0)) - sum(lgamma(priors$alpha0)) +
    sum((priors$alpha0 - 1) * Eq_log_pi)

#print_named(Eq_log_prior_pi, "res-v2.txt")

  # Expected log p(mu | Lambda) and Expected log p(Lambda).

  Eq_log_prior_mu_given_Lambda <- 0
  Eq_log_prior_Lambda <- K * log_wishart_constant(priors$W0, priors$nu0,
    hDnu0log2, qDDm1logpi, sum_hlg_n0p1mseqD)

  for (k in seq_len(K))
  {
    Mk <- D / pars$mu_betas[k] +
      pars$Lambda_nus[[k]] * quadratic_byrows(pars$mu_ms[k,,drop=FALSE], priors$m0, pars$Lambda_Ws[[k]])

    Eq_log_prior_mu_given_Lambda <- Eq_log_prior_mu_given_Lambda +
      hDlogbeta0 - hDlog2pi + 0.5 * Eq_log_Lambda[k] - hbeta0 * Mk

    Eq_log_prior_Lambda <- Eq_log_prior_Lambda +
      hnu0mdp1 * Eq_log_Lambda[k] - hnus[k] * sum(priors$W0_inv * pars$Lambda_Ws[[k]])

    # The equivalence sum(A * t(B)) equals sum(diag(A %*% B))
    # is used in sum(priors$W0_inv * t(pars$Lambda_Ws[[k]])).
    # Because of symmetry of the matrices the transpose can be avoided.
    # The first option is much faster.
  }

#print_named(Eq_log_prior_mu_given_Lambda, "res-v2.txt")
#print_named(Eq_log_prior_Lambda, "res-v2.txt")

  # Entropy terms.

  entropy_q_Z <- -sum(rprobs * log(pmax(rprobs, .Machine$double.eps)))
  entropy_q_pi <- entropy_dirichlet(pars$pi_alphas, K)

  entropy_q_mu_Lambda <- 0

  for (k in seq_len(K))
  {
    entropy_mu_given_Lambda <- hD_1plog2pi  -
      hD * log(pars$mu_betas[k]) - 0.5 * Eq_log_Lambda[k]

    entropy_q_mu_Lambda <- entropy_q_mu_Lambda + entropy_mu_given_Lambda +
      entropy_wishart(pars$Lambda_Ws[[k]], pars$Lambda_nus[k],
                      D, hD, Dlog2, hDnulog2[k], qDDm1logpi, sum_hlg_np1mseqD[k], hnumDm1[k])
  }

  elbo <- Eq_log_lik + Eq_log_Z + Eq_log_prior_pi +
    Eq_log_prior_mu_given_Lambda + Eq_log_prior_Lambda +
    entropy_q_Z + entropy_q_pi + entropy_q_mu_Lambda

#print_named(entropy_q_Z, "res-v2.txt")
#print_named(entropy_q_pi, "res-v2.txt")
#print_named(entropy_mu_given_Lambda, "res-v2.txt")
#print_named(entropy_q_mu_Lambda, "res-v2.txt")
#print_named(elbo, "res-v2.txt")

  if (debug) return(
    list(Eq_log_lik = Eq_log_lik, Eq_log_Z = Eq_log_Z,
         Eq_log_prior_pi = Eq_log_prior_pi,
         Eq_log_prior_mu_given_Lambda = Eq_log_prior_mu_given_Lambda,
         Eq_log_prior_Lambda = Eq_log_prior_Lambda,
         entropy_q_Z = entropy_q_Z, entropy_q_pi = entropy_q_pi, entropy_q_mu_Lambda = entropy_q_mu_Lambda,
         elbo = elbo))

  return(elbo)
}


## -------------------------------------------------------------------
## Posterior summary. Estimated means.
## -------------------------------------------------------------------

mog_nw_posteriors <- function(model)
{
  pars <- model$pars
  K <- model$K
  D <- model$D

  # Mixture probabilities: E_q[pi_k].
  pi_hat <- pars$pi_alphas / sum(pars$pi_alphas)

  # Component means: E_q[mu_k].
  # mu_hat <- pars$mu_ms

  # Precision and covariance summaries.
  Lambda_hat <- vector("list", K)
  Sigma_hat_plugin <- vector("list", K)
  Sigma_hat_mean <- vector("list", K)
  Sigma_mu_mean <- vector("list", K)

  for (k in seq_len(K))
  {
    # E_q[Lambda_k].
    Lambda_hat[[k]] <- pars$Lambda_nus[k] * pars$Lambda_Ws[[k]]

    # Plug-in covariance from the expected precision:
    # solve(E_q[Lambda_k]).
    Sigma_hat_plugin[[k]] <- solve(Lambda_hat[[k]])

    # Mean covariance E_q[Lambda_k^{-1}], if it exists.
    #
    # If Lambda_k ~ Wishart(W_k, nu_k), then
    # E[Lambda_k^{-1}] = W_k^{-1} / (nu_k - D - 1),
    # provided nu_k > D + 1.
    if (pars$Lambda_nus[k] > D + 1)
    {
      Sigma_hat_mean[[k]] <- solve(pars$Lambda_Ws[[k]]) / (pars$Lambda_nus[k] - D - 1)

      # Since q(mu_k | Lambda_k) has covariance
      # (beta_k Lambda_k)^(-1), the posterior mean covariance
      # of mu_k conditional uncertainty is:
      Sigma_mu_mean[[k]] <- Sigma_hat_mean[[k]] / pars$mu_betas[k]
    } else {
      Sigma_hat_mean[[k]] <- matrix(NA_real_, D, D)
      Sigma_mu_mean[[k]] <- matrix(NA_real_, D, D)
    }
  }

  list(pis = pi_hat, mus = pars$mu_ms,
    Lambda = Lambda_hat, Sigma_plugin = Sigma_hat_plugin,
    Sigma_mean = Sigma_hat_mean, Sigma_mu_mean = Sigma_mu_mean)
}


## -------------------------------------------------------------------
## Print method for object of class "mog_cavi".
## -------------------------------------------------------------------

print.mog_normal_wishart_cavi <- function(x, digits = 4, ...)
{
  if (!inherits(x, "mog_normal_wishart_cavi")) {
    stop("'x' must be an object of class 'mog_normal_wishart_cavi'.")
  }

  model <- x$model

  if (is.null(x$posteriors)) {
    posteriors <- mog_nw_posteriors(model)
  } else {
    posteriors <- x$posteriors
  }

  cat("\n")
  cat("Variational Bayes fit for mixture-of-Gaussians model\n")
  cat("====================================================\n")

  cat("\nModel dimensions\n")
  cat("----------------------------------------------------------------\n")
  cat("Number of observations:", model$N, "\n")
  cat("Number of variables D:", model$D, "\n")
  cat("Number of mixture components K:", model$K, "\n")

  cat("\nConvergence\n")
  cat("----------------------------------------------------------------\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$niter, "\n")

  if (!is.null(x$elbo_path))
  {
    elbo_used <- x$elbo_path[is.finite(x$elbo_path)]

    if (length(elbo_used) > 0)
    {
      cat("Initial ELBO:", round(elbo_used[1], digits), "\n")
      cat("Final ELBO:", round(tail(elbo_used, 1), digits), "\n")
    }
  }

  cat("\nPosteriors summary (means)\n")
  cat("----------------------------------------------------------------\n")

  cat("\nMixture probabilities:\n")
  print(round(posteriors$pis, digits))

  cat("\nComponent means, m_k:\n")
  print(round(posteriors$mus, digits))

  cat("\nComponent plug-in covariances solve(E[Lambda_k])\n")
  invisible(lapply(posteriors$Sigma_plugin, function(z) { print(round(z, digits)) }))

  cat("\nComponent mean covariances E[Lambda_k^{-1}]\n")
  invisible(lapply(posteriors$Sigma_mean, function(z) {
    print(round(z, digits))
  }))

  cat("\nPosterior mean covariance of component means E[Lambda_k^{-1}] / beta_k\n")
  invisible(lapply(posteriors$Sigma_mu_mean, function(z) {
    print(round(z, digits))
  }))

  cat("\n")
  invisible(x)
}


## -------------------------------------------------------------------
## Normal-Wishart family.
## -------------------------------------------------------------------

# Using update* naming convention would be neater, e.g.:
# update_auxiliary, update, update_responsibilities.
# However, using names not starting by the same substring can
# make list searches by name slightly easier/faster.

mog_normal_wishart <- list(
  update = mog_nw_upd,
  resp_update = mog_nw_upd_resp,
  init = mog_nw_init,
  elbo = mog_nw_elbo,
  posteriors = mog_nw_posteriors,
  name = "mog_normal_wishart"
)
