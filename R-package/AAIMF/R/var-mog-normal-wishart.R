## Variational Bayes for the
## VAR mixture-of-Gaussians-noise model.


# -------------------------------------------------------------------
# Initialization of parameters.
# -------------------------------------------------------------------

var_mog_nw_init <- function(model,
  init_method = c("residual_kmeans", "random"),
  kmeans_nstart = 10, kmeans_prob = 0.95,
  ridge_B = 1e-6, ridge_cov = 1e-6,
  seed = NULL,
  update_resp = FALSE)
{
  # Build the model with make_model() for preliminary checks.

  init_method <- match.arg(init_method)

  X <- model$X
  Xlags <- model$Xlags
  priors <- model$priors
  N <- model$N
  K <- model$K
  D <- model$D
  Dp <- ncol(Xlags)
  M <- ncol(Xlags) * D

  # ------------------------------------------------------------
  # Fixed prior parameters.
  # ------------------------------------------------------------

  if (length(priors$m0) == 1L)
    priors$m0 <- rep(priors$m0, D)

  if (length(priors$m0) != D)
    stop("'priors$m0' must have length D.")

  if (is.null(priors$W0))
    priors$W0 <- diag(D)

  if (!is.matrix(priors$W0))
    priors$W0 <- as.matrix(priors$W0)

  if (!all(dim(priors$W0) == c(D, D)))
    stop("'priors$W0' must be a D x D matrix.")

  if (is.null(priors$nu0))
    priors$nu0 <- D + 2

  if (priors$nu0 <= D - 1)
    stop("'priors$nu0' must be larger than D - 1.")

  if (priors$kappa0 <= 0)
    stop("'priors$kappa0' must be positive.")

  if (priors$alpha_a0 <= 0 || priors$alpha_r0 <= 0)
    stop("'alpha_a0' and 'alpha_r0' must be positive.")

  if (priors$lambda0 <= 0)
    stop("'priors$lambda0' must be positive.")

  if (ridge_B < 0 || ridge_cov < 0)
    stop("'ridge_B' and 'ridge_cov' must be non-negative.")

  # Not needed. This is done by make_model() which ensure proper arguments
  #priors <- modifyList.topLevel(formals(family$init)[[1]], priors)
  #default_priors <- formals(var_mog_nw_init)[[1]]
  #priors <- modifyList_priorsPars(formals(var_mog_nw_init)[[1]], priors)

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

  #if (K == 1L)
  #{
  #  rprobs <- matrix(1.0, nrow = N, ncol = 1L)
  #} else
  if (init_method == "residual_kmeans")
  {
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
  # Initialize Normal-Wishart parameters.
  # ------------------------------------------------------------

  m <- matrix(NA_real_, K, D)
  kappa <- rep(NA_real_, K)
  W <- vector("list", K)
  nu <- rep(NA_real_, K)

  W0_inv <- solve(priors$W0)

  for (k in seq_len(K))
  {
    gamma_k <- rprobs[,k]

    ebar_k <- as.vector(crossprod(gamma_k, E_init) / Nk_safe[k])

    S_k <- matrix(0.0, D, D)

    for (n in seq_len(N))
    {
      diff_nk <- E_init[n, ] - ebar_k
      S_k <- S_k + gamma_k[n] * tcrossprod(diff_nk)
    }

    kappa[k] <- priors$kappa0 + Nk[k]
    nu[k] <- priors$nu0 + Nk[k]

    m[k,] <- (priors$kappa0 * priors$m0 + Nk[k] * ebar_k) / (priors$kappa0 + Nk[k])

    diff_prior <- ebar_k - priors$m0

    W_inv_k <- W0_inv + S_k +
      (priors$kappa0 * Nk[k] / (priors$kappa0 + Nk[k])) * tcrossprod(diff_prior)

    W_inv_k <- W_inv_k + ridge_cov * diag(D)

    W[[k]] <- solve(W_inv_k)

  } # end loop seq_len(K)

  colnames(m) <- colnames(X)

  # ------------------------------------------------------------
  # Initialize q(alpha).
  # ------------------------------------------------------------

  # Start with a simple diagonal covariance for q(B).
  Sigma_B_initial <- diag(M)

  bbar <- as.vector(Bbar)
  alpha_a <- priors$alpha_a0 + M / 2

  alpha_r <- priors$alpha_r0 + 0.5 * (sum(bbar^2) + sum(diag(Sigma_B_initial)))

  Eq_alpha <- alpha_a / alpha_r

  # ------------------------------------------------------------
  # Refine initial Sigma_B using the current expected precisions.
  # ------------------------------------------------------------

  Prec_B <- Eq_alpha * diag(M)

  for (k in seq_len(K))
  {
    Lambda_bar_k <- nu[k] * W[[k]]

    H_weighted <- sqrt(rprobs[, k]) * Xlags
    S_hh_k <- crossprod(H_weighted)

    Prec_B <- Prec_B + kronecker(Lambda_bar_k, S_hh_k)
  }

  Prec_B <- Prec_B + ridge_B * diag(M)
  Sigma_B <- solve(Prec_B)

  # Recompute alpha_r with the refined Sigma_B.
  alpha_r <- priors$alpha_r0 + 0.5 * (sum(bbar^2) + sum(diag(Sigma_B)))

  # ------------------------------------------------------------
  # Collect parameters in a list.
  # ------------------------------------------------------------

  pars <- list(lambdas = lambdas, Bbar = Bbar, bbar = bbar,
    Sigma_B = Sigma_B, Prec_B = Prec_B,
    alpha_a = alpha_a, alpha_r = alpha_r,
    m = m, kappa = kappa,
    W = W, nu = nu, Nk = Nk)

  # (Optional) Update responsibilities using the current
  # variational parameters.

  if (isTRUE(update_resp))
  {
    rprobs <- var_mog_nw_upd_resp(model, debug = FALSE)

    pars$Nk <- colSums(probs)
    pars$lambdas <- priors$lambda0 + pars$Nk
  }

  # In R, lists and most objects are passed by value with copy-on-modify.
  # R typically avoids copying the untouched elements due to to copy-on-write semantics.
  # So it is efficient to return the whole list (with updated and the unmodified elements).

  #list(priors = priors, rprobs = rprobs, pars = pars)
  model$priors <- priors
  model$rprobs <- rprobs
  model$pars <- pars

  model
}


# -------------------------------------------------------------------
# CAVI updates.
# -------------------------------------------------------------------

var_mog_make_G <- function(h, D)
{
  # h is a vector of length Dp.
  # Returns G such that B^T h = G vec(B),
  # where vec(B) follows R's column-major order.
  kronecker(diag(D), matrix(h, nrow = 1L))
}

# quadratic_byrows <- function(X, A)
# {
#   # Computes x_i^T A x_i for every row x_i of X.
#   rowSums((X %*% A) * X)
# }

var_mog_nw_upd <- function(model)
{
  # Variational update for the Generalized VAR model with
  # mixture-of-Gaussian innovations and Normal-Wishart priors.
  #
  # X      : N x D response matrix, already lag-trimmed.
  # Xlags  : N x Dp lag matrix.
  # rprobs : N x K responsibility matrix gamma_nk.
  #
  # priors should contain:
  # lambda0, alpha_a0, alpha_r0,
  # m0, kappa0, W0, nu0.
  #
  # pars should contain current:
  # Bbar, Sigma_B, alpha_a, alpha_r,
  # m, kappa, W, nu.

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

  #if (ncol(X) != D)
  #  stop("'D' is not equal to ncol(X).")
  #
  #if (nrow(Xlags) != N)
  #  stop("'X' and 'Xlags' must have the same number of rows.")
  #
  #if (nrow(rprobs) != N || ncol(rprobs) != K)
  #  stop("'rprobs' must be an N x K matrix.")

  m0 <- priors$m0
  if (length(m0) == 1L)
    m0 <- rep(m0, D)

  W0 <- priors$W0
  W0_inv <- solve(W0)

  # Current expectations.

  Eq_alpha <- pars$alpha_a / pars$alpha_r

  Lambda_bar <- vector("list", K)
  for (k in seq_len(K))
    Lambda_bar[[k]] <- pars$nu[k] * pars$W[[k]]

  # Effective component sizes.

  colSums_rprobs <- colSums(rprobs)
  Nk_safe <- pmax(colSums_rprobs, .Machine$double.eps)

  # Update q(pi).

  lambdas <- priors$lambda0 + colSums_rprobs

  # Update q(B) = N(vec(Bbar), Sigma_B).

  Prec_B <- Eq_alpha * diag(M)
  lin_B_mat <- matrix(0.0, Dp, D)

  for (k in seq_len(K))
  {
    Lambda_k <- Lambda_bar[[k]]

    # Sum_n gamma_nk h_n h_n^T.
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

  # Update q(alpha) = Gamma(alpha_a, alpha_r).

  alpha_a <- priors$alpha_a0 + M / 2

  alpha_r <- priors$alpha_r0 + 0.5 * (sum(bbar^2) + sum(diag(Sigma_B)))

  # Update q(mu_k, Lambda_k), Normal-Wishart.

  Ebar <- X - Xlags %*% Bbar

  # C_n = G_n Sigma_B G_n^T accounts for posterior uncertainty in B.
  C_list <- vector("list", N)

  for (n in seq_len(N))
  {
    G_n <- var_mog_make_G(Xlags[n, ], D)
    C_list[[n]] <- G_n %*% Sigma_B %*% t(G_n)
  }

  m <- matrix(NA_real_, nrow = K, ncol = D)
  kappa <- rep(NA_real_, K)
  W <- vector("list", K)
  nu <- rep(NA_real_, K)

  for (k in seq_len(K))
  {
    gamma_k <- rprobs[,k]

    if (colSums_rprobs[k] <= .Machine$double.eps)
    {
      ebar_k <- m0
    } else {
      ebar_k <- as.vector(crossprod(gamma_k, Ebar) / Nk_safe[k])
    }

    S_k <- matrix(0.0, nrow = D, ncol = D)

    for (n in seq_len(N))
    {
      diff_nk <- Ebar[n, ] - ebar_k

      S_k <- S_k +
        gamma_k[n] *
        (tcrossprod(diff_nk) + C_list[[n]])
    }

    kappa[k] <- priors$kappa0 + colSums_rprobs[k]
    nu[k] <- priors$nu0 + colSums_rprobs[k]

    m[k, ] <- (priors$kappa0 * m0 + colSums_rprobs[k] * ebar_k) / (priors$kappa0 + colSums_rprobs[k])

    diff_prior <- ebar_k - m0

    W_inv_k <- W0_inv + S_k +
      (priors$kappa0 * colSums_rprobs[k]) / (priors$kappa0 + colSums_rprobs[k]) * tcrossprod(diff_prior)

    W[[k]] <- solve(W_inv_k)
  }

  list(lambdas = lambdas,
    Bbar = Bbar, bbar = bbar, Prec_B = Prec_B,
    Sigma_B = Sigma_B,
    alpha_a = alpha_a, alpha_r = alpha_r,
    m = m, kappa = kappa,
    W = W, nu = nu)
}

var_mog_nw_upd_resp <- function(model, debug = TRUE)
{
  # Variational responsibility update for the VAR-MoG model.

  X <- model$X
  Xlags <- model$Xlags
  pars <- model$pars
  K <- model$K
  N <- model$N
  D <- model$D

  Eq_log_pi <- digamma(pars$lambdas) - digamma(sum(pars$lambdas))

  Eq_log_Lambda <- rep(NA_real_, K)
  Lambda_bar <- vector("list", K)

  for (k in seq_len(K))
  {
    logdet_W_k <- as.numeric(determinant(pars$W[[k]], logarithm = TRUE)$modulus)

    Eq_log_Lambda[k] <- sum(digamma((pars$nu[k] + 1L - seq_len(D)) / 2)) +
      D * log(2) + logdet_W_k

    Lambda_bar[[k]] <- pars$nu[k] * pars$W[[k]]
  }

  Xhat <- Xlags %*% pars$Bbar

  # C_n = G_n Sigma_B G_n^T.
  C_list <- vector("list", N)

  for (n in seq_len(N))
  {
    G_n <- var_mog_make_G(Xlags[n, ], D)
    C_list[[n]] <- G_n %*% pars$Sigma_B %*% t(G_n)
  }

  Q <- matrix(NA_real_, nrow = N, ncol = K)

  for (k in seq_len(K))
  {
    Delta_k <- sweep(X - Xhat, 2L, pars$m[k, ], "-",
      check.margin = FALSE)

    #quad_k <- quadratic_byrows(Delta_k, Lambda_bar[[k]])
    # Computes x_i^T A x_i for every row x_i of X.
    quad_k <- rowSums((Delta_k %*% Lambda_bar[[k]]) * Delta_k)

    trace_k <- rep(NA_real_, N)
    for (n in seq_len(N))
    {
      trace_k[n] <- sum(diag(Lambda_bar[[k]] %*% C_list[[n]]))
    }

    Q[, k] <- quad_k + trace_k + D / pars$kappa[k]
  }

  scores <- matrix(Eq_log_pi, nrow = N, ncol = K, byrow = TRUE) +
    0.5 * matrix(Eq_log_Lambda, nrow = N, ncol = K, byrow = TRUE) -
    0.5 * D * log(2 * pi) - 0.5 * Q

  rprobs <- softmax_rows(scores)

  if (debug && !isTRUE(all.equal(rowSums(rprobs), rep(1, N))))
    stop("Responsibilities do not sum to one.")

  rprobs
}


# -------------------------------------------------------------------
# Evidence lower bound (ELBO) for the MoG model.
# Analytical ELBO.
# -------------------------------------------------------------------

# FIXME compare with as.numeric(determinant(A, logarithm = TRUE)A$modulus)
logdet_spd_mat_v2 <- function(A)
{
  R <- chol(A)
  2 * sum(log(diag(R)))
}

log_multivariate_gamma <- function(a, D)
{
  # log Gamma_D(a)
  D * (D - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(D)) / 2))
}

var_mog_nw_elbo <- function(model,
  #X,
  #priors, rprobs, pars, etas,
  #D, K, Xlags,
  debug = TRUE)
{
  # ELBO for the Generalized VAR model with mixture-of-Gaussian
  # innovations and Normal-Wishart priors.

  X <- model$X
  Xlags <- model$Xlags
  priors <- model$priors
  rprobs <- model$rprobs
  pars <- model$pars
  #etas <- model$etas
  K <- model$K
  D <- model$D

  N <- nrow(X)
  Dp <- ncol(Xlags)
  M <- Dp * D

  log2pi <- log(2 * pi)

  m0 <- priors$m0
  if (length(m0) == 1L)
    m0 <- rep(m0, D)

  W0 <- priors$W0
  W0_inv <- solve(W0)
  logdet_W0 <- logdet_spd_mat_v2(W0)

  # ------------------------------------------------------------------
  # Auxiliary elements.
  # ------------------------------------------------------------------

  # FIXME define once pars$xxx instead of searching pars$xxx every time

  Lambda_sum <- sum(pars$lambdas)
  digamma_Lambda_sum <- digamma(Lambda_sum)
  digamma_lambdas <- digamma(pars$lambdas)

  Eq_log_pi <- digamma_lambdas - digamma_Lambda_sum

  Eq_log_alpha <- digamma(pars$alpha_a) - log(pars$alpha_r)
  Eq_alpha <- pars$alpha_a / pars$alpha_r

  bbar <- as.vector(pars$Bbar)

  Eq_btb <- sum(bbar^2) + sum(diag(pars$Sigma_B))

  Eq_log_Lambda <- rep(NA_real_, K)
  Eq_Lambda <- vector("list", K)

  for (k in seq_len(K))
  {
    logdet_W_k <- logdet_spd_mat_v2(pars$W[[k]])

    Eq_log_Lambda[k] <- sum(digamma((pars$nu[k] + 1 - seq_len(D)) / 2)) +
      D * log(2) + logdet_W_k

    Eq_Lambda[[k]] <- pars$nu[k] * pars$W[[k]]
  }

  # ------------------------------------------------------------------
  # Expected likelihood term.
  # ------------------------------------------------------------------

  Xhat <- Xlags %*% pars$Bbar

  C_list <- vector("list", N)

  for (n in seq_len(N))
  {
    G_n <- var_mog_make_G(Xlags[n, ], D)
    C_list[[n]] <- G_n %*% pars$Sigma_B %*% t(G_n)
  }

  Qnk <- matrix(NA_real_, nrow = N, ncol = K)

  for (k in seq_len(K))
  {
    Delta_k <- sweep(X - Xhat, 2L, pars$m[k, ], "-",
      check.margin = FALSE)

    #quad_k <- quadratic_byrows(Delta_k, Eq_Lambda[[k]])
    # Computes x_i^T A x_i for every row x_i of X.
    quad_k <- rowSums((Delta_k %*% Eq_Lambda[[k]]) * Delta_k)

    trace_k <- rep(NA_real_, N)

    for (n in seq_len(N))
    {
      trace_k[n] <- sum(diag(Eq_Lambda[[k]] %*% C_list[[n]]))
    }

    Qnk[, k] <- quad_k + trace_k + D / pars$kappa[k]
  }

  Eq_log_lik <- rprobs * (
    matrix(Eq_log_pi, nrow = N, ncol = K, byrow = TRUE) +
    0.5 * matrix(Eq_log_Lambda, nrow = N, ncol = K, byrow = TRUE) -
    0.5 * D * log2pi - 0.5 * Qnk)

  # ------------------------------------------------------------------
  # Expected log-prior terms.
  # ------------------------------------------------------------------

  # Dirichlet prior for pi.
  Eq_log_prior_pi <- lgamma(K * priors$lambda0) -
    K * lgamma(priors$lambda0) + (priors$lambda0 - 1) * sum(Eq_log_pi)

  # Gaussian prior for vec(B) | alpha.
  Eq_log_prior_B <- 0.5 * M * Eq_log_alpha - 0.5 * M * log2pi -
    0.5 * Eq_alpha * Eq_btb

  # Gamma prior for alpha.
  Eq_log_prior_alpha <- priors$alpha_a0 * log(priors$alpha_r0) -
    lgamma(priors$alpha_a0) + (priors$alpha_a0 - 1) * Eq_log_alpha -
    priors$alpha_r0 * Eq_alpha

  # Normal prior for mu_k | Lambda_k.
  Eq_log_prior_mu <- 0.0

  for (k in seq_len(K))
  {
    diff_mu_k <- pars$m[k, ] - m0

    Q_mu0_k <- as.numeric(t(diff_mu_k) %*% Eq_Lambda[[k]] %*% diff_mu_k) +
      D / pars$kappa[k]

    Eq_log_prior_mu <- Eq_log_prior_mu +
      0.5 * D * log(priors$kappa0) + 0.5 * Eq_log_Lambda[k] -
      0.5 * D * log2pi - 0.5 * priors$kappa0 * Q_mu0_k
  }

  # Wishart prior for Lambda_k.
  Eq_log_prior_Lambda <- 0.0

  for (k in seq_len(K))
  {
    Eq_log_prior_Lambda <- Eq_log_prior_Lambda +
      0.5 * (priors$nu0 - D - 1) * Eq_log_Lambda[k] -
      0.5 * sum(diag(W0_inv %*% Eq_Lambda[[k]])) -
      0.5 * priors$nu0 * D * log(2) -
      0.5 * priors$nu0 * logdet_W0 -
      log_multivariate_gamma(priors$nu0 / 2, D)
  }

  # ------------------------------------------------------------------
  # Entropy terms: -E_q[log q].
  # ------------------------------------------------------------------

  # Entropy of q(pi), Dirichlet.
  entropy_q_pi <- sum(lgamma(pars$lambdas)) -
    lgamma(Lambda_sum) + (Lambda_sum - K) * digamma_Lambda_sum -
    sum((pars$lambdas - 1) * digamma_lambdas)

  # Entropy of q(B), multivariate normal.
  logdet_Sigma_B <- logdet_spd_mat_v2(pars$Sigma_B)

  entropy_q_B <- 0.5 * M * (log2pi + 1) + 0.5 * logdet_Sigma_B

  # Entropy of q(alpha), Gamma(shape, rate).
  entropy_q_alpha <- pars$alpha_a - log(pars$alpha_r) +
    lgamma(pars$alpha_a) + (1 - pars$alpha_a) * digamma(pars$alpha_a)

  # Entropy of q(mu_k | Lambda_k).
  # This is -E_q[log q(mu_k | Lambda_k)].
  entropy_q_mu_given_Lambda <- 0.0

  for (k in seq_len(K))
  {
    entropy_q_mu_given_Lambda <- entropy_q_mu_given_Lambda +
      0.5 * D * (log2pi + 1) - 0.5 * D * log(pars$kappa[k]) -
      0.5 * Eq_log_Lambda[k]
  }

  # Entropy of q(Lambda_k), Wishart.
  entropy_q_Lambda <- 0.0

  for (k in seq_len(K))
  {
    logdet_W_k <- logdet_spd_mat_v2(pars$W[[k]])

    # Directly compute -E_q[log q(Lambda_k)].
    Eq_log_q_Lambda_k <- 0.5 * (pars$nu[k] - D - 1) * Eq_log_Lambda[k] -
      0.5 * sum(diag(solve(pars$W[[k]]) %*% Eq_Lambda[[k]])) -
      0.5 * pars$nu[k] * D * log(2) -
      0.5 * pars$nu[k] * logdet_W_k -
      log_multivariate_gamma(pars$nu[k] / 2, D)

    entropy_q_Lambda <- entropy_q_Lambda - Eq_log_q_Lambda_k
  }

  # Entropy of q(Z): categorical allocations.
  entropy_q_Z <- -sum(rprobs * log(pmax(rprobs, .Machine$double.eps)))

  # ------------------------------------------------------------------
  # ELBO.
  # ------------------------------------------------------------------

  elbo <- sum(Eq_log_lik) + Eq_log_prior_pi + Eq_log_prior_B +
    Eq_log_prior_alpha + Eq_log_prior_mu + Eq_log_prior_Lambda +
    entropy_q_pi + entropy_q_B + entropy_q_alpha +
    entropy_q_mu_given_Lambda + entropy_q_Lambda + entropy_q_Z

  if (debug)
  {
    return(list(Eq_log_lik = Eq_log_lik,
        Qnk = Qnk,
        Eq_log_pi = Eq_log_pi, Eq_log_alpha = Eq_log_alpha,
        Eq_alpha = Eq_alpha,
        Eq_log_Lambda = Eq_log_Lambda, Eq_Lambda = Eq_Lambda,
        Eq_btb = Eq_btb,
        Eq_log_prior_pi = Eq_log_prior_pi, Eq_log_prior_B = Eq_log_prior_B,
        Eq_log_prior_alpha = Eq_log_prior_alpha, Eq_log_prior_mu = Eq_log_prior_mu,
        Eq_log_prior_Lambda = Eq_log_prior_Lambda,
        entropy_q_pi = entropy_q_pi, entropy_q_B = entropy_q_B,
        entropy_q_alpha = entropy_q_alpha,
        entropy_q_mu_given_Lambda = entropy_q_mu_given_Lambda,
        entropy_q_Lambda = entropy_q_Lambda, entropy_q_Z = entropy_q_Z,
        elbo = elbo))
  }

  elbo
}


# -------------------------------------------------------------------
# Posterior summary. Estimated means.
# -------------------------------------------------------------------

B_to_A <- function(B, D, p)
{
  A <- vector("list", p)

  for (i in seq_len(p))
  {
    rows_i <- ((i - 1L) * D + 1L):(i * D)
    A[[i]] <- t(B[rows_i,,drop = FALSE])
  }

  A
}

var_mog_nw_posteriors <- function(model)
{
  # FIXME see reuse mog_nw_posteriors()

  pars <- model$pars
  K <- model$K
  D <- model$D
  p <- model$p

  # mixture probabilities
  pi_hat <- pars$lambdas / sum(pars$lambdas)

  # innovation means
  # mu_hat <- pars$m

  Lambda_hat <- vector("list", K)
  Sigma_hat_plugin <- vector("list", K)
  Sigma_hat_mean <- vector("list", K)

  for (k in seq_len(K))
  {
    # E_q[Lambda_k]
    Lambda_hat[[k]] <- pars$nu[k] * pars$W[[k]]

    # Innovation's plug-in covariance from the expected precision.
    Sigma_hat_plugin[[k]] <- solve(Lambda_hat[[k]])

    # Innovation's mean covariances E[Lambda_k^{-1} (if it exists).
    # If Lambda ~ Wishart(W, nu), then:
    # E[Lambda^{-1}] = W^{-1} / (nu - D - 1).
    if (pars$nu[k] > D + 1)
    {
      Sigma_hat_mean[[k]] <- solve(pars$W[[k]]) / (pars$nu[k] - D - 1)
    } else {
      Sigma_hat_mean[[k]] <- matrix(NA, D, D)
    }
  }

  # VAR coefficient matrices (Bbar converted to A matrices)
  A_hat <- B_to_A(pars$Bbar, D, p)

  list(pi = pi_hat,
       mu = pars$m, Lambda = Lambda_hat,
       Sigma_plugin = Sigma_hat_plugin, Sigma_mean = Sigma_hat_mean,
       A = A_hat)
}


# -------------------------------------------------------------------
# Print method for object of class "var_mog_normal_wishart_cavi".
# -------------------------------------------------------------------

print.var_mog_normal_wishart_cavi <- function(x, digits = 4, ...)
{
  if (!inherits(x, "var_mog_normal_wishart_cavi")) {
    stop("'x' must be an object of class 'var_mog_normal_wishart_cavi'.")
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
    posteriors <- var_mog_nw_posteriors(model)
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

  cat("\nInnovation plug-in covariances solve(E[Lambda_k])\n")
  invisible(lapply(posteriors$Sigma_plugin, function(x) print(round(x, digits))))
  cat("\nInnovation mean covariances E[Lambda_k^{-1}]\n")
  invisible(lapply(posteriors$Sigma_mean, function(x) print(round(x, digits))))

  cat("\nVAR coefficient matrices (Bbar converted to A matrices)\n")
  invisible(lapply(posteriors$A, function(x) print(round(x, digits))))

  cat("\n")
  invisible(x)
}


# -------------------------------------------------------------------
# Vector Autoregressive Normal-Wishart family.
# -------------------------------------------------------------------

var_mog_normal_wishart <- list(
  update = var_mog_nw_upd,
  resp_update = var_mog_nw_upd_resp,
  init = var_mog_nw_init,
  elbo = var_mog_nw_elbo,
  posteriors = var_mog_nw_posteriors,
  name = "var_mog_normal_wishart"
)
