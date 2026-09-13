## Variational Bayes for the VAR mixture-of-Gaussians-noise model
## with diagonal innovation covariance matrices and Normal-Gamma priors.
##
## Model:
##   X[n, ] = B^T h_n + epsilon_n,
##   epsilon_n | z_n = k ~ N_D(mu_k, diag(lambda_k)^(-1)).
##
## Variational family:
##   q(Z) q(pi) q(vec(B)) q(alpha)
##   prod_k prod_d q(mu_kd, lambda_kd),
##
## where
##   q(mu_kd | lambda_kd) = N(m_kd, (kappa_k lambda_kd)^(-1)),
##   q(lambda_kd)         = Gamma(a_kd, b_kd),
##
## and all Gamma distributions use the shape-rate parameterization.


# -------------------------------------------------------------------
# Utilities.
# -------------------------------------------------------------------

var_mog_ng_expand_D <- function(x, D, name)
{
  if (length(x) == 1L)
    x <- rep(x, D)

  if (length(x) != D)
    stop(sprintf("'%s' must have length 1 or D.", name))

  as.numeric(x)
}

var_mog_ng_logdet_spd <- function(A)
{
  R <- chol(A)
  2 * sum(log(diag(R)))
}

var_mog_ng_make_G <- function(h, D)
{
  # If b = vec(B), using R's column-major vectorization, then
  # B^T h = G(h) b.
  kronecker(diag(D), matrix(h, nrow = 1L))
}


# -------------------------------------------------------------------
# Initialization of parameters.
# -------------------------------------------------------------------

var_mog_ng_init <- function(model,
  init_method = c("residual_kmeans", "random"),
  kmeans_nstart = 10, kmeans_prob = 0.95,
  ridge_B = 1e-6, ridge_var = 1e-8,
  seed = NULL,
  update_resp = FALSE)
{
  # Build the model with make_model() before calling this function so
  # dimensions and the common model fields are already available.

  init_method <- match.arg(init_method)

  X <- model$X
  Xlags <- model$Xlags
  priors <- model$priors
  N <- model$N
  K <- model$K
  D <- model$D
  Dp <- ncol(Xlags)
  M <- Dp * D

  # ------------------------------------------------------------
  # Fixed prior parameters.
  # ------------------------------------------------------------

  priors$m0 <- var_mog_ng_expand_D(priors$m0, D, "priors$m0")
  priors$a0 <- var_mog_ng_expand_D(priors$a0, D, "priors$a0")
  priors$b0 <- var_mog_ng_expand_D(priors$b0, D, "priors$b0")

  if (priors$kappa0 <= 0)
    stop("'priors$kappa0' must be positive.")

  if (any(priors$a0 <= 0) || any(priors$b0 <= 0))
    stop("'priors$a0' and 'priors$b0' must be positive.")

  if (priors$alpha_a0 <= 0 || priors$alpha_r0 <= 0)
    stop("'alpha_a0' and 'alpha_r0' must be positive.")

  if (priors$lambda0 <= 0)
    stop("'priors$lambda0' must be positive.")

  if (ridge_B < 0 || ridge_var < 0)
    stop("'ridge_B' and 'ridge_var' must be non-negative.")

  if (K < 1L)
    stop("'K' must be at least one.")

  if (K == 1L && init_method == "residual_kmeans")
    kmeans_prob <- 1

  if (K > 1L && (kmeans_prob <= 1 / K || kmeans_prob >= 1))
    stop("For K > 1, 'kmeans_prob' must lie in (1/K, 1).")

  # ------------------------------------------------------------
  # Initial ridge VAR fit.
  # ------------------------------------------------------------

  XtX <- crossprod(Xlags)
  XtY <- crossprod(Xlags, X)

  Bbar <- solve(XtX + ridge_B * diag(Dp), XtY)
  E_init <- X - Xlags %*% Bbar

  # ------------------------------------------------------------
  # Initial responsibilities.
  # ------------------------------------------------------------

  if (K == 1L)
  {
    rprobs <- matrix(1, N, 1L)

  } else if (init_method == "residual_kmeans") {
    km <- kmeans(E_init, centers = K, nstart = kmeans_nstart)
    rprobs <- matrix((1 - kmeans_prob) / (K - 1), N, K)
    rprobs[cbind(seq_len(N), km$cluster)] <- kmeans_prob

  } else {
    if (!is.null(seed))
      set.seed(seed)

    rprobs <- matrix(rgamma(N * K, shape = 1, rate = 1), N, K)
    rprobs <- sweep(rprobs, 1L, rowSums(rprobs), "/")
  }

  # ------------------------------------------------------------
  # Effective component sizes and q(pi).
  # ------------------------------------------------------------

  Nk <- colSums(rprobs)
  Nk_safe <- pmax(Nk, .Machine$double.eps)
  lambdas <- priors$lambda0 + Nk

  # ------------------------------------------------------------
  # Initialize Normal-Gamma parameters.
  # ------------------------------------------------------------

  m <- matrix(NA_real_, K, D)
  kappa <- priors$kappa0 + Nk
  a <- matrix(NA_real_, K, D)
  b <- matrix(NA_real_, K, D)

  for (k in seq_len(K))
  {
    gamma_k <- rprobs[, k]

    if (Nk[k] <= .Machine$double.eps)
    {
      ebar_k <- priors$m0
    } else {
      ebar_k <- as.vector(crossprod(gamma_k, E_init) / Nk_safe[k])
    }

    m[k, ] <- (priors$kappa0 * priors$m0 + Nk[k] * ebar_k) /
      (priors$kappa0 + Nk[k])

    a[k, ] <- priors$a0 + 0.5 * Nk[k]

    for (d in seq_len(D))
    {
      scatter_kd <- sum(gamma_k * (E_init[, d] - ebar_k[d])^2)
      prior_shift <- priors$kappa0 * Nk[k] /
        (priors$kappa0 + Nk[k]) * (ebar_k[d] - priors$m0[d])^2

      b[k, d] <- priors$b0[d] + 0.5 * (scatter_kd + prior_shift) + ridge_var
    }
  }

  colnames(m) <- colnames(X)
  colnames(a) <- colnames(X)
  colnames(b) <- colnames(X)

  # ------------------------------------------------------------
  # Initialize q(alpha).
  # ------------------------------------------------------------

  Sigma_B_initial <- diag(M)
  bbar <- as.vector(Bbar)

  alpha_a <- priors$alpha_a0 + M / 2
  alpha_r <- priors$alpha_r0 +
    0.5 * (sum(bbar^2) + sum(diag(Sigma_B_initial)))

  Eq_alpha <- alpha_a / alpha_r

  # ------------------------------------------------------------
  # Refine initial Sigma_B using current expected precisions.
  # ------------------------------------------------------------

  Prec_B <- Eq_alpha * diag(M)

  for (k in seq_len(K))
  {
    lambda_bar_k <- a[k, ] / b[k, ]
    H_weighted <- sqrt(rprobs[, k]) * Xlags
    S_hh_k <- crossprod(H_weighted)

    Prec_B <- Prec_B + kronecker(diag(lambda_bar_k), S_hh_k)
  }

  Prec_B <- Prec_B + ridge_B * diag(M)
  Sigma_B <- solve(Prec_B)

  alpha_r <- priors$alpha_r0 +
    0.5 * (sum(bbar^2) + sum(diag(Sigma_B)))

  # ------------------------------------------------------------
  # Collect parameters.
  # ------------------------------------------------------------

  pars <- list( lambdas = lambdas,
    Bbar = Bbar, bbar = bbar,
    Sigma_B = Sigma_B, Prec_B = Prec_B,
    alpha_a = alpha_a, alpha_r = alpha_r,
    m = m,
    kappa = kappa,
    a = a, b = b,
    Nk = Nk)

  model$priors <- priors
  model$rprobs <- rprobs
  model$pars <- pars

  # Optional responsibility refinement must use a model containing
  # the newly initialized variational parameters.
  if (isTRUE(update_resp))
  {
    model$rprobs <- var_mog_ng_upd_resp(model, debug = FALSE)
    model$pars$Nk <- colSums(model$rprobs)
    model$pars$lambdas <- priors$lambda0 + model$pars$Nk
  }

  model
}


# -------------------------------------------------------------------
# CAVI parameter updates.
# -------------------------------------------------------------------

var_mog_ng_upd <- function(model)
{
  # One coordinate-ascent update for the generalized VAR model with
  # mixture-of-Gaussian innovations and Normal-Gamma priors.

  X <- model$X
  Xlags <- model$Xlags
  priors <- model$priors
  rprobs <- model$rprobs
  pars <- model$pars
  K <- model$K
  N <- model$N
  D <- model$D

  Dp <- ncol(Xlags)
  M <- Dp * D

  m0 <- var_mog_ng_expand_D(priors$m0, D, "priors$m0")
  a0 <- var_mog_ng_expand_D(priors$a0, D, "priors$a0")
  b0 <- var_mog_ng_expand_D(priors$b0, D, "priors$b0")

  # Current expectations.

  Eq_alpha <- pars$alpha_a / pars$alpha_r
  lambda_bar <- pars$a / pars$b

  # Effective component sizes.

  Nk <- colSums(rprobs)
  Nk_safe <- pmax(Nk, .Machine$double.eps)

  # ------------------------------------------------------------
  # Update q(pi).
  # ------------------------------------------------------------

  lambdas <- priors$lambda0 + Nk

  # ------------------------------------------------------------
  # Update q(B) = N(vec(Bbar), Sigma_B).
  # ------------------------------------------------------------

  Prec_B <- Eq_alpha * diag(M)
  lin_B_mat <- matrix(0.0, Dp, D)

  for (k in seq_len(K))
  {
    Lambda_k <- diag(lambda_bar[k, ], D, D)

    H_weighted <- sqrt(rprobs[, k]) * Xlags
    S_hh_k <- crossprod(H_weighted)

    Prec_B <- Prec_B + kronecker(Lambda_k, S_hh_k)

    X_minus_mk <- sweep(X, 2L, pars$m[k, ], "-", check.margin = FALSE)
    X_weighted <- sweep(X_minus_mk, 1L, rprobs[, k], "*",
      check.margin = FALSE)

    lin_B_mat <- lin_B_mat + crossprod(Xlags, X_weighted) %*% Lambda_k
  }

  Sigma_B <- solve(Prec_B)
  bbar <- as.vector(Sigma_B %*% as.vector(lin_B_mat))
  Bbar <- matrix(bbar, Dp, D)

  # ------------------------------------------------------------
  # Update q(alpha) = Gamma(alpha_a, alpha_r).
  # ------------------------------------------------------------

  alpha_a <- priors$alpha_a0 + M / 2
  alpha_r <- priors$alpha_r0 +
    0.5 * (sum(bbar^2) + sum(diag(Sigma_B)))

  # ------------------------------------------------------------
  # Update q(mu_kd, lambda_kd), Normal-Gamma.
  # ------------------------------------------------------------

  Ebar <- X - Xlags %*% Bbar

  # C_n = G_n Sigma_B G_n^T accounts for posterior uncertainty in B.
  # Only diag(C_n) enters a diagonal-precision Normal-Gamma update.
  C_diag <- matrix(NA_real_, N, D)

  for (n in seq_len(N))
  {
    G_n <- var_mog_ng_make_G(Xlags[n, ], D)
    C_n <- G_n %*% Sigma_B %*% t(G_n)
    C_diag[n, ] <- diag(C_n)
  }

  m <- matrix(NA_real_, K, D)
  kappa <- priors$kappa0 + Nk
  a <- matrix(NA_real_, K, D)
  b <- matrix(NA_real_, K, D)

  for (k in seq_len(K))
  {
    gamma_k <- rprobs[, k]

    if (Nk[k] <= .Machine$double.eps)
    {
      ebar_k <- m0
    } else {
      ebar_k <- as.vector(crossprod(gamma_k, Ebar) / Nk_safe[k])
    }

    m[k, ] <- (priors$kappa0 * m0 + Nk[k] * ebar_k) /
      (priors$kappa0 + Nk[k])

    a[k, ] <- a0 + 0.5 * Nk[k]

    for (d in seq_len(D))
    {
      # Expected within-component residual scatter.  The C_diag term
      # propagates uncertainty in B into the precision update.
      scatter_kd <- sum(gamma_k *
        ((Ebar[, d] - ebar_k[d])^2 + C_diag[, d]))

      prior_shift <- priors$kappa0 * Nk[k] /
        (priors$kappa0 + Nk[k]) * (ebar_k[d] - m0[d])^2

      b[k, d] <- b0[d] + 0.5 * (scatter_kd + prior_shift)
    }
  }

  colnames(m) <- colnames(X)
  colnames(a) <- colnames(X)
  colnames(b) <- colnames(X)

  list(
    lambdas = lambdas,
    Bbar = Bbar,
    bbar = bbar,
    Prec_B = Prec_B,
    Sigma_B = Sigma_B,
    alpha_a = alpha_a,
    alpha_r = alpha_r,
    m = m,
    kappa = kappa,
    a = a,
    b = b,
    Nk = Nk
  )
}

# -------------------------------------------------------------------
# Responsibility update.
# -------------------------------------------------------------------

var_mog_ng_upd_resp <- function(model, debug = TRUE)
{
  X <- model$X
  Xlags <- model$Xlags
  pars <- model$pars
  K <- model$K
  N <- model$N
  D <- model$D

  Eq_log_pi <- digamma(pars$lambdas) - digamma(sum(pars$lambdas))
  Eq_log_lambda <- digamma(pars$a) - log(pars$b)
  Eq_lambda <- pars$a / pars$b

  Xhat <- Xlags %*% pars$Bbar

  # Diagonal entries of C_n = G_n Sigma_B G_n^T.
  C_diag <- matrix(NA_real_, N, D)

  for (n in seq_len(N))
  {
    G_n <- var_mog_ng_make_G(Xlags[n, ], D)
    C_n <- G_n %*% pars$Sigma_B %*% t(G_n)
    C_diag[n, ] <- diag(C_n)
  }

  Q <- matrix(NA_real_, N, K)

  for (k in seq_len(K))
  {
    Delta_k <- sweep(X - Xhat, 2L, pars$m[k, ], "-",
      check.margin = FALSE)

    # E[sum_d lambda_kd (x_nd - B^T h_n - mu_kd)^2].
    Q[, k] <- rowSums(
      sweep(Delta_k^2 + C_diag, 2L, Eq_lambda[k, ], "*",
        check.margin = FALSE)
    ) + D / pars$kappa[k]
  }

  scores <- matrix(Eq_log_pi, N, K, byrow = TRUE) +
    0.5 * matrix(rowSums(Eq_log_lambda), N, K, byrow = TRUE) -
    0.5 * D * log(2 * pi) -
    0.5 * Q

  rprobs <- softmax_rows(scores)

  if (debug && !isTRUE(all.equal(rowSums(rprobs), rep(1, N))))
    stop("Responsibilities do not sum to one.")

  rprobs
}

# -------------------------------------------------------------------
# Analytical evidence lower bound.
# -------------------------------------------------------------------

var_mog_ng_elbo <- function(model, debug = TRUE)
{
  X <- model$X
  Xlags <- model$Xlags
  priors <- model$priors
  rprobs <- model$rprobs
  pars <- model$pars
  K <- model$K
  D <- model$D

  N <- nrow(X)
  Dp <- ncol(Xlags)
  M <- Dp * D

  log2pi <- log(2 * pi)

  m0 <- var_mog_ng_expand_D(priors$m0, D, "priors$m0")
  a0 <- var_mog_ng_expand_D(priors$a0, D, "priors$a0")
  b0 <- var_mog_ng_expand_D(priors$b0, D, "priors$b0")

  # ------------------------------------------------------------
  # Useful expectations.
  # ------------------------------------------------------------

  lambda_sum <- sum(pars$lambdas)
  digamma_lambda_sum <- digamma(lambda_sum)
  digamma_lambdas <- digamma(pars$lambdas)

  Eq_log_pi <- digamma_lambdas - digamma_lambda_sum

  Eq_log_alpha <- digamma(pars$alpha_a) - log(pars$alpha_r)
  Eq_alpha <- pars$alpha_a / pars$alpha_r

  bbar <- as.vector(pars$Bbar)
  Eq_btb <- sum(bbar^2) + sum(diag(pars$Sigma_B))

  Eq_log_lambda <- digamma(pars$a) - log(pars$b)
  Eq_lambda <- pars$a / pars$b

  # ------------------------------------------------------------
  # Expected likelihood term.
  # ------------------------------------------------------------

  Xhat <- Xlags %*% pars$Bbar

  C_diag <- matrix(NA_real_, N, D)

  for (n in seq_len(N))
  {
    G_n <- var_mog_ng_make_G(Xlags[n, ], D)
    C_n <- G_n %*% pars$Sigma_B %*% t(G_n)
    C_diag[n, ] <- diag(C_n)
  }

  Qnk <- matrix(NA_real_, N, K)

  for (k in seq_len(K))
  {
    Delta_k <- sweep(X - Xhat, 2L, pars$m[k, ], "-",
      check.margin = FALSE)

    Qnk[, k] <- rowSums(
      sweep(Delta_k^2 + C_diag, 2L, Eq_lambda[k, ], "*",
        check.margin = FALSE)
    ) + D / pars$kappa[k]
  }

  Eq_log_lik <- rprobs * (
    matrix(Eq_log_pi, N, K, byrow = TRUE) +
    0.5 * matrix(rowSums(Eq_log_lambda), N, K, byrow = TRUE) -
    0.5 * D * log2pi -
    0.5 * Qnk
  )

  # ------------------------------------------------------------
  # Expected log-prior terms.
  # ------------------------------------------------------------

  # Symmetric Dirichlet prior for pi.
  Eq_log_prior_pi <- lgamma(K * priors$lambda0) -
    K * lgamma(priors$lambda0) +
    (priors$lambda0 - 1) * sum(Eq_log_pi)

  # Gaussian prior for vec(B) | alpha.
  Eq_log_prior_B <- 0.5 * M * Eq_log_alpha -
    0.5 * M * log2pi -
    0.5 * Eq_alpha * Eq_btb

  # Gamma prior for alpha.
  Eq_log_prior_alpha <- priors$alpha_a0 * log(priors$alpha_r0) -
    lgamma(priors$alpha_a0) +
    (priors$alpha_a0 - 1) * Eq_log_alpha -
    priors$alpha_r0 * Eq_alpha

  # Normal prior for mu_kd | lambda_kd.
  Eq_log_prior_mu <- 0.0

  for (k in seq_len(K))
  {
    Q_mu0_k <- Eq_lambda[k, ] * (pars$m[k, ] - m0)^2 +
      1 / pars$kappa[k]

    Eq_log_prior_mu <- Eq_log_prior_mu +
      0.5 * D * log(priors$kappa0) -
      0.5 * D * log2pi +
      0.5 * sum(Eq_log_lambda[k, ]) -
      0.5 * priors$kappa0 * sum(Q_mu0_k)
  }

  # Gamma priors for coordinate precisions lambda_kd.
  Eq_log_prior_lambda <- 0.0

  for (k in seq_len(K))
  {
    Eq_log_prior_lambda <- Eq_log_prior_lambda + sum(
      a0 * log(b0) - lgamma(a0) +
      (a0 - 1) * Eq_log_lambda[k, ] -
      b0 * Eq_lambda[k, ]
    )
  }

  # ------------------------------------------------------------
  # Entropy terms: -E_q[log q].
  # ------------------------------------------------------------

  entropy_q_pi <- sum(lgamma(pars$lambdas)) -
    lgamma(lambda_sum) +
    (lambda_sum - K) * digamma_lambda_sum -
    sum((pars$lambdas - 1) * digamma_lambdas)

  logdet_Sigma_B <- var_mog_ng_logdet_spd(pars$Sigma_B)
  entropy_q_B <- 0.5 * M * (log2pi + 1) + 0.5 * logdet_Sigma_B

  entropy_q_alpha <- pars$alpha_a - log(pars$alpha_r) +
    lgamma(pars$alpha_a) +
    (1 - pars$alpha_a) * digamma(pars$alpha_a)

  # Conditional Gaussian entropy for q(mu_kd | lambda_kd).
  entropy_q_mu_given_lambda <- 0.0

  for (k in seq_len(K))
  {
    entropy_q_mu_given_lambda <- entropy_q_mu_given_lambda +
      0.5 * D * (log2pi + 1) -
      0.5 * D * log(pars$kappa[k]) -
      0.5 * sum(Eq_log_lambda[k, ])
  }

  # Gamma entropy for q(lambda_kd), shape-rate form.
  entropy_q_lambda <- sum(
    pars$a - log(pars$b) + lgamma(pars$a) +
    (1 - pars$a) * digamma(pars$a)
  )

  entropy_q_Z <- -sum(rprobs * log(pmax(rprobs, .Machine$double.eps)))

  # ------------------------------------------------------------
  # ELBO.
  # ------------------------------------------------------------

  elbo <- sum(Eq_log_lik) +
    Eq_log_prior_pi + Eq_log_prior_B + Eq_log_prior_alpha +
    Eq_log_prior_mu + Eq_log_prior_lambda +
    entropy_q_pi + entropy_q_B + entropy_q_alpha +
    entropy_q_mu_given_lambda + entropy_q_lambda + entropy_q_Z

  if (debug)
  {
    return(list(
      Eq_log_lik = Eq_log_lik,
      Qnk = Qnk,
      Eq_log_pi = Eq_log_pi,
      Eq_log_alpha = Eq_log_alpha,
      Eq_alpha = Eq_alpha,
      Eq_log_lambda = Eq_log_lambda,
      Eq_lambda = Eq_lambda,
      Eq_btb = Eq_btb,
      Eq_log_prior_pi = Eq_log_prior_pi,
      Eq_log_prior_B = Eq_log_prior_B,
      Eq_log_prior_alpha = Eq_log_prior_alpha,
      Eq_log_prior_mu = Eq_log_prior_mu,
      Eq_log_prior_lambda = Eq_log_prior_lambda,
      entropy_q_pi = entropy_q_pi,
      entropy_q_B = entropy_q_B,
      entropy_q_alpha = entropy_q_alpha,
      entropy_q_mu_given_lambda = entropy_q_mu_given_lambda,
      entropy_q_lambda = entropy_q_lambda,
      entropy_q_Z = entropy_q_Z,
      elbo = elbo
    ))
  }

  elbo
}

# -------------------------------------------------------------------
# Posterior summaries.
# -------------------------------------------------------------------

var_mog_ng_posteriors <- function(model)
{
  pars <- model$pars
  K <- model$K
  D <- model$D
  p <- model$p

  pi_hat <- pars$lambdas / sum(pars$lambdas)

  lambda_hat <- pars$a / pars$b
  Sigma_plugin <- 1 / lambda_hat

  # For lambda ~ Gamma(a, b), E[lambda^{-1}] = b / (a - 1),
  # provided a > 1.
  Sigma_mean <- matrix(NA_real_, K, D)
  ok <- pars$a > 1
  Sigma_mean[ok] <- pars$b[ok] / (pars$a[ok] - 1)

  A_hat <- B_to_A(pars$Bbar, D, p)

  list(pi = pi_hat,
      mu = pars$m,
      lambda = lambda_hat,
      Sigma_plugin_diag = Sigma_plugin, Sigma_mean_diag = Sigma_mean,
      A = A_hat,
      alpha = pars$alpha_a / pars$alpha_r)
}

# -------------------------------------------------------------------
# Print method for object of class "var_mog_normal_gamma_cavi".
# -------------------------------------------------------------------

print.var_mog_normal_gamma_cavi <- function(x, digits = 4, ...)
{
  # TODO report alpha

  if (!inherits(x, "var_mog_normal_gamma_cavi")) {
    stop("'x' must be an object of class 'var_mog_normal_gamma_cavi'.")
  }

  fmt <- function(z) {
    if (is.numeric(z)) {
      round(z, digits = digits)
    } else {
      z
    }
  }

  model <- x$model

  if (is.null(x$posteriors)) {
    posteriors <- var_mog_ng_posteriors(model)
  } else {
    posteriors <- x$posteriors
  }

  cat("\n")
  cat("Variational Bayes fit for VAR mixture-of-Gaussians-noise model\n")
  cat("==============================================================\n")

  cat("\nModel dimensions\n")
  cat("----------------------------------------------------------------\n")
  cat("Number of observations used:", model$N, "\n")
  cat("AR order p:", model$p, "\n")
  cat("Number of mixture components K:", model$K, "\n")

  cat("\nConvergence\n")
  cat("----------------------------------------------------------------\n")

  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$niter, "\n")
  if (!is.null(x$elbo_path))
  {
    cat("Initial ELBO:", round(x$elbo_path[1], digits), "\n")
    elbo_used <- x$elbo_path[is.finite(x$elbo_path)]
    if (length(elbo_used) > 0) {
      cat("Final ELBO:", round(tail(elbo_used, 1), digits), "\n")
    }
  }

  cat("\nPosteriors summary (means)\n")
  cat("----------------------------------------------------------------\n")

  cat("\nMixture probabilities:", round(posteriors$pi, digits), "\n")

  cat("\nInnovation means, m_k:\n")
  print(round(posteriors$mu, digits))

  cat("\nComponent expected precisions E[lambda_kd]:\n")
  print(round(posteriors$lambda, digits))

  cat("\nInnovation plug-in variances 1 / E[lambda_kd]\n")
  invisible(lapply(posteriors$Sigma_plugin_diag, function(x) print(round(x, digits))))
  cat("\nInnovation mean variances E[1 / lambda_kd]\n")
  invisible(lapply(posteriors$Sigma_mean_diag, function(x) print(round(x, digits))))

  cat("\nVAR coefficient matrices (Bbar converted to A matrices)\n")
  invisible(lapply(posteriors$A, function(x) print(round(x, digits))))

  cat("\n")
  invisible(x)
}


# -------------------------------------------------------------------
# Vector autoregressive Normal-Gamma family.
# -------------------------------------------------------------------

var_mog_normal_gamma <- list(
  update = var_mog_ng_upd,
  resp_update = var_mog_ng_upd_resp,
  init = var_mog_ng_init,
  elbo = var_mog_ng_elbo,
  posteriors = var_mog_ng_posteriors,
  name = "var_mog_normal_gamma"
)
