## Variational Bayes for the
## Normal-Gamma mixture-of-Gaussians-noise model.


## -------------------------------------------------------------------
## Initialization of parameters.
## -------------------------------------------------------------------

mog_ng_init <- function(model)
{
  # Build the model with make_model() for preliminary checks.

  # In this model, only priors need to be initialized here.
  # The remaining parameters are initialized by the general-purpose function mog_init().

  priors <- model$priors
  K <- model$K
  D <- model$D

  if (length(priors$alpha0) == 1) {
    priors$alpha0 <- rep(priors$alpha0, K)
  } else
  if (length(priors$alpha0) != K)
    stop("'priors$alpha0' must have length 1 or K.")

  if (length(priors$beta0) != 1 || priors$beta0 <= 0)
    stop("'priors$beta0' must be a positive scalar.")

  if (is.null(priors$m0)) {
    priors$m0 <- rep(0, D)
  } else
  if (length(priors$m0) != D)
    stop("'priors$m0' must have length D.")

  if (length(priors$a0) == 1) {
    priors$a0 <- rep(priors$a0, D)
  } else
  if (length(priors$a0) != D)
    stop("'priors$a0' must have length 1 or D.")

  if (length(priors$b0) == 1) {
    priors$b0 <- rep(priors$b0, D)
  } else
  if (length(priors$b0) != D)
    stop("'priors$b0' must have length 1 or D.")

  if (any(priors$alpha0 <= 0)) stop("'priors$alpha0' must be positive.")
  if (any(priors$a0 <= 0)) stop("'priors$a0' must be positive.")
  if (any(priors$b0 <= 0)) stop("'priors$b0' must be positive.")

  # Not needed. This is done by make_model() which ensure proper arguments
  #priors <- modifyList_priorsPars(formals(mog_ng_init)[[1]], priors)

  #list(fixed_pars = fixed_pars, pars = NULL)
  #list(priors = priors, rprobs = NULL, pars = NULL) #init_info = list()
  model$priors <- priors
  model
}


## -------------------------------------------------------------------
## CAVI updates.
## -------------------------------------------------------------------

mog_ng_upd_aux <- function(X, rprobs, colSums_rprobs, K, D)
{
  # Internal function. Assumes callers pass valid arguments.

  # Unlike in the Normal-Wishart, here full D x D matrix is
  # not needed, only diagonal responsibility-weighted terms.

  eps <- 1e-12

  xbar <- matrix(0, K, D)
  S <- matrix(0, K, D)

  for (k in seq_len(K))
  {
    if (colSums_rprobs[k] > eps)
    {
      xbar[k,] <- colSums(rprobs[,k] * X) / colSums_rprobs[k]
      #Xcentered <- sweep(X, MARGIN = 2, STATS = xbar[k,], FUN = "-", check.margin = FALSE)
      Xcentered <- sweep(X, 2, xbar[k,], "-", FALSE)
      S[k,] <- colSums(rprobs[,k] * Xcentered^2) / colSums_rprobs[k]
    } # else (initialized with zeros)
  }

  list(xbar = xbar, S = S)
}

mog_ng_upd <- function(model)
{
  # Internal function. Assumes callers pass valid arguments.

  X <- model$X
  priors <- model$priors
  rprobs <- model$rprobs
  K <- model$K
  D <- model$D

  colSums_rprobs <- colSums(rprobs)
  aux <- mog_ng_upd_aux(X, rprobs, colSums_rprobs, K, D)

  pi_alphas <- priors$alpha0 + colSums_rprobs

  mu_betas <- priors$beta0 + colSums_rprobs
  mu_ms <- (priors$beta0 * matrix(priors$m0, K, D, byrow = TRUE) +
    matrix(colSums_rprobs, K, D, byrow = FALSE) * aux$xbar) / mu_betas

  weights <- priors$beta0 * colSums_rprobs / mu_betas

  lambda_as <- matrix(priors$a0, K, D, byrow = TRUE) +
    matrix(0.5 * colSums_rprobs, K, D, byrow = FALSE)

  lambda_bs <- matrix(priors$b0, K, D, byrow = TRUE) +
    0.5 * matrix(colSums_rprobs, K, D, byrow = FALSE) * aux$S +
    0.5 * matrix(weights, K, D, byrow = FALSE) *
      sweep(aux$xbar, MARGIN = 2L, STATS = priors$m0, FUN = "-", check.margin = FALSE)^2

  list(pi_alphas = pi_alphas, mu_ms = mu_ms, mu_betas = mu_betas,
       lambda_as = lambda_as, lambda_bs = lambda_bs)
}

mog_ng_upd_resp <- function(model, debug = TRUE)
{
  # Internal function. Assumes callers pass valid arguments.
  # Recall: Amortization is not addressed by CAVI.

  X <- model$X
  pars <- model$pars
  K <- model$K
  N <- model$N
  D <- model$D

  Eq_log_pi <- digamma(pars$pi_alphas) - digamma(sum(pars$pi_alphas))
  Eq_log_lambda <- digamma(pars$lambda_as) - log(pars$lambda_bs)
  Eq_lambda <- pars$lambda_as / pars$lambda_bs

  hELmDlog2pi <- 0.5 * rowSums(Eq_log_lambda) - 0.5 * D * log(2 * pi)

  scores <- matrix(0, N, K)

  for (k in seq_len(K))
  {
    #Xcentered <- sweep(X, MARGIN = 2L, STATS = pars$mu_ms[k,], FUN = "-", check.margin = FALSE)
    Xcentered <- sweep(X, 2L, pars$mu_ms[k,], "-", FALSE)
    Qnk <- D / pars$mu_betas[k] +
      rowSums(sweep(Xcentered^2, 2L, Eq_lambda[k,], "*", FALSE))
      #rowSums(sweep(Xcentered^2, MARGIN = 2L, STATS = Eq_lambda[k,], FUN = "*", check.margin = FALSE))

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

entropy_gamma_rate <- function(a, b)
{
  # Entropy of Gamma(a,b), shape-rate parameterization.
  a - log(b) + lgamma(a) + (1 - a) * digamma(a)
}

mog_ng_elbo <- function(model, debug = TRUE)
{
  X <- model$X
  priors <- model$priors
  rprobs <- model$rprobs
  pars <- model$pars
  K <- model$K
  D <- model$D

  # Constants and expectations.

  log2pi <- log(2 * pi)
  hDlog2pi <- 0.5 * D * log2pi
  hlogbeta0 <- 0.5 * log(priors$beta0)
  hlog2pi <- 0.5 * log2pi
  hbeta0 <- 0.5 * priors$beta0

  Eq_log_pi <- digamma(pars$pi_alphas) - digamma(sum(pars$pi_alphas))
  Eq_log_lambda <- digamma(pars$lambda_as) - log(pars$lambda_bs)
  Eq_lambda <- pars$lambda_as / pars$lambda_bs

  # Expected log likelihood E_q[log p(X | Z, mu, lambda)].

  Eq_log_lik <- 0
  for (k in seq_len(K))
  {
    #Xcentered <- sweep(X, MARGIN = 2L, STATS = pars$mu_ms[k,], FUN = "-", check.margin = FALSE)
    Xcentered <- sweep(X, 2L, pars$mu_ms[k,], "-", FALSE)
    Qnk <- D / pars$mu_betas[k] +
      rowSums(sweep(Xcentered^2, 2L, Eq_lambda[k,], "*", FALSE))
      #rowSums(sweep(Xcentered^2, MARGIN = 2L, STATS = Eq_lambda[k,], FUN = "*", check.margin = FALSE))

    Eq_log_lik <- Eq_log_lik +
      sum(rprobs[,k] * (0.5 * sum(Eq_log_lambda[k,]) - hDlog2pi - 0.5 * Qnk))
  }

  # Expected log p(Z | pi).

  Eq_log_Z <- sum(sweep(x = rprobs, MARGIN = 2L,
    STATS = Eq_log_pi, FUN = "*", check.margin = FALSE))

  # Expected log p(pi).

  Eq_log_prior_pi <- lgamma(sum(priors$alpha0)) - sum(lgamma(priors$alpha0)) +
    sum((priors$alpha0 - 1) * Eq_log_pi)

  # Expected log p(mu | lambda) and Expected log p(lambda).

  a0_mat <- matrix(priors$a0, K, D, byrow = TRUE)
  b0_mat <- matrix(priors$b0, K, D, byrow = TRUE)
  m0_mat <- matrix(priors$m0, K, D, byrow = TRUE)

  Mkd <- 1 / matrix(pars$mu_betas, K, D, byrow = FALSE) +
    Eq_lambda * (pars$mu_ms - m0_mat)^2

  Eq_log_prior_mu_given_lambda <- sum(hlogbeta0 - hlog2pi +
    0.5 * Eq_log_lambda - hbeta0 * Mkd)

  Eq_log_prior_lambda <- sum(a0_mat * log(b0_mat) - lgamma(a0_mat) +
    (a0_mat - 1) * Eq_log_lambda - b0_mat * Eq_lambda)

  # Entropy terms.

  entropy_q_Z <- -sum(rprobs * log(pmax(rprobs, .Machine$double.eps)))
  entropy_q_pi <- entropy_dirichlet(pars$pi_alphas, K)

  entropy_mu_given_lambda <- 0.5 * (1 + log2pi) -
    0.5 * log(matrix(pars$mu_betas, K, D, byrow = FALSE)) -
    0.5 * Eq_log_lambda

  entropy_q_mu_lambda <- sum(entropy_mu_given_lambda +
    entropy_gamma_rate(pars$lambda_as, pars$lambda_bs))

  elbo <- Eq_log_lik + Eq_log_Z + Eq_log_prior_pi +
    Eq_log_prior_mu_given_lambda + Eq_log_prior_lambda +
    entropy_q_Z + entropy_q_pi + entropy_q_mu_lambda

  if (debug) return(
    list(Eq_log_lik = Eq_log_lik, Eq_log_Z = Eq_log_Z,
         Eq_log_prior_pi = Eq_log_prior_pi,
         Eq_log_prior_mu_given_lambda = Eq_log_prior_mu_given_lambda,
         Eq_log_prior_lambda = Eq_log_prior_lambda,
         entropy_q_Z = entropy_q_Z, entropy_q_pi = entropy_q_pi,
         entropy_q_mu_lambda = entropy_q_mu_lambda,
         elbo = elbo))

  return(elbo)
}


## -------------------------------------------------------------------
## Posterior summary. Estimated means.
## -------------------------------------------------------------------

mog_ng_posteriors <- function(model)
{
  if (is.null(model$pars))
    stop("'model' must contain 'pars'.")

  pars <- model$pars
  K <- model$K
  D <- model$D

  # Mixture probabilities: E_q[pi_k].
  pi_hat <- pars$pi_alphas / sum(pars$pi_alphas)

  # Component means: E_q[mu_k].
  mu_hat <- pars$mu_ms

  # Component precisions: E_q[lambda_kd].
  #
  # In the Normal-Gamma model each coordinate has its own precision,
  # so lambda_as and lambda_bs are K x D matrices.
  lambda_hat <- pars$lambda_as / pars$lambda_bs

  # Plug-in variances from the expected precision:
  # 1 / E_q[lambda_kd].
  sigma2_plugin <- 1 / lambda_hat

  # Mean variances E_q[1 / lambda_kd], if they exist.
  #
  # If lambda_kd ~ Gamma(a_kd, b_kd), shape-rate, then:
  # E[1 / lambda_kd] = b_kd / (a_kd - 1), for a_kd > 1.
  sigma2_mean <- matrix(NA_real_, nrow = K, ncol = D)

  valid_inv_moment <- pars$lambda_as > 1
  sigma2_mean[valid_inv_moment] <- pars$lambda_bs[valid_inv_moment] /
    (pars$lambda_as[valid_inv_moment] - 1)

  # Posterior conditional variance of component means.
  #
  # q(mu_kd | lambda_kd) = N(mu_ms[k,d], (mu_betas[k] lambda_kd)^(-1)).
  #
  # Therefore:
  # E[Var(mu_kd | lambda_kd)] = E[1 / lambda_kd] / mu_betas[k].
  mu_var_mean <- sigma2_mean /
    matrix(pars$mu_betas, nrow = K, ncol = D, byrow = FALSE)

  list(pis = pi_hat, mus = mu_hat,
    #Lambdas = unname(as.list(as.data.frame(lambda_hat))),
    # NOTE review for NormalWishart: one element per mixture component (not one list element per matrix colum)
    Lambdas = lapply(seq_len(K), function(k) lambda_hat[k,]),
    # TODO keep because used in previous versions (see confint);
    # in the future stick to report Lambdas for easier comparison with EM output
    lambda = lambda_hat,
    sigma2_plugin = sigma2_plugin, sigma2_mean = sigma2_mean,
    mu_var_mean = mu_var_mean)
}


## -------------------------------------------------------------------
## Print method for object of class "mog_normal_gamma_cavi".
## -------------------------------------------------------------------

print.mog_normal_gamma_cavi <- function(x, digits = 4, ...)
{
  if (!inherits(x, "mog_normal_gamma_cavi")) {
    stop("'x' must be an object of class 'mog_normal_gamma_cavi'.")
  }

  model <- x$model

  if (is.null(x$posteriors)) {
    posteriors <- mog_ng_posteriors(model)
  } else {
    posteriors <- x$posteriors
  }

  cat("\n")
  cat("Variational Bayes fit for mixture-of-Gaussians model\n")
  cat("Normal-Gamma prior with diagonal component precisions\n")
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

  cat("\nComponent expected precisions E[lambda_kd]:\n")
  print(round(posteriors$lambda, digits))

  cat("\nComponent plug-in variances 1 / E[lambda_kd]:\n")
  print(round(posteriors$sigma2_plugin, digits))

  cat("\nComponent mean variances E[1 / lambda_kd]:\n")
  print(round(posteriors$sigma2_mean, digits))

  cat("\nPosterior mean variances of component means E[1 / lambda_kd] / beta_k:\n")
  print(round(posteriors$mu_var_mean, digits))

  cat("\n")
  invisible(x)
}


## -------------------------------------------------------------------
## Confidence intervals for class mog_normal_gamma_cavi
## -------------------------------------------------------------------

confint.mog_normal_gamma_cavi <- function(object, level = 0.95, struct = NULL)
{
  if (!inherits(object, "mog_normal_gamma_cavi")) {
    stop("'object' must be an object of class 'mog_normal_gamma_cavi'.")
  }

  stopifnot(level >= 0 && level <= 1)
  #stop("'level' must be a number between 0 and 1.")

  model <- object$model
  pars <- model$pars
  pi_alphas <- pars$pi_alphas
  mu_ms <- pars$mu_ms
  mu_betas <- pars$mu_betas
  lambda_as <- pars$lambda_as
  lambda_bs <- pars$lambda_bs
  K <- model$K
  D <- model$D

  lower_prob <- (1 - level) / 2
  upper_prob <- (1 + level) / 2

  # Posterior point summaries.
  post <- mog_ng_posteriors(model)

  # ------------------------------------------------------------
  # Mixture probabilities pi_k.
  #
  # q(pi) = Dirichlet(pi_alphas).
  # Marginally:
  # pi_k ~ Beta(pi_alphas[k], sum(pi_alphas) - pi_alphas[k]).
  # ------------------------------------------------------------

  alpha_hat <- sum(pi_alphas)

  pi_lower <- qbeta(lower_prob, shape1 = pi_alphas, shape2 = alpha_hat - pi_alphas)
  pi_upper <- qbeta(upper_prob, shape1 = pi_alphas, shape2 = alpha_hat - pi_alphas)

  pi_ci <- cbind(lower = pi_lower, estimate = post$pis, upper = pi_upper)
  rownames(pi_ci) <- paste0("k", seq_len(K))

  npars <- length(pi_lower)

  # ------------------------------------------------------------
  # Component means mu_kd.
  #
  # q(mu_kd | lambda_kd) = N(m_kd, (beta_k lambda_kd)^(-1))
  # q(lambda_kd) = Gamma(a_kd, b_kd), shape-rate.
  #
  # Marginally:
  # mu_kd ~ t_{2 a_kd}
  #   location = m_kd
  #   scale    = sqrt(b_kd / (a_kd beta_k)).
  # ------------------------------------------------------------

  mu_beta_mat <- matrix(mu_betas, nrow = K, ncol = D, byrow = FALSE)
  mu_scale <- sqrt(lambda_bs / (lambda_as * mu_beta_mat))
  mu_lower <- mu_ms + qt(lower_prob, df = 2 * lambda_as) * mu_scale
  mu_upper <- mu_ms + qt(upper_prob, df = 2 * lambda_as) * mu_scale

  npars <- npars + prod(dim(mu_lower))

  # ------------------------------------------------------------
  # Component precisions lambda_kd.
  #
  # q(lambda_kd) = Gamma(a_kd, b_kd), shape-rate.
  # ------------------------------------------------------------

  lambda_lower <- qgamma(lower_prob, shape = lambda_as, rate = lambda_bs)
  lambda_upper <- qgamma(upper_prob, shape = lambda_as, rate = lambda_bs)

  npars <- npars + prod(dim(lambda_lower))

  # ------------------------------------------------------------
  # Component variances sigma2_kd = 1 / lambda_kd.
  #
  # Equal-tailed intervals are obtained by transforming Gamma
  # quantiles. Since sigma2 = 1 / lambda is decreasing in lambda,
  # the lower and upper quantiles are reversed.
  # ------------------------------------------------------------

  sigma2_lower <- 1 / qgamma(upper_prob, shape = lambda_as, rate = lambda_bs)
  sigma2_upper <- 1 / qgamma(lower_prob, shape = lambda_as, rate = lambda_bs)

  # ------------------------------------------------------------
  # Component standard deviations sigma_kd.
  # ------------------------------------------------------------

  sigma_lower <- sqrt(sigma2_lower)
  sigma_upper <- sqrt(sigma2_upper)
  sigma_est_plugin <- sqrt(post$sigma2_plugin)
  sigma_est_mean <- sqrt(post$sigma2_mean)

  # output

  #res <- list(
  #  level = level,
  #  pis = pi_ci,
  #
  #  mus = list(
  #    estimate = post$mus,
  #    lower = mu_lower,
  #    upper = mu_upper,
  #    scale = mu_scale,
  #    df = 2 * lambda_as,
  #    variance_mean = post$mu_var_mean),
  #
  #  Lambdas = list(
  #    estimate = post$lambda,
  #    lower = lambda_lower,
  #    upper = lambda_upper),
  #
  #  sigma2 = list(
  #    estimate_plugin = post$sigma2_plugin,
  #    estimate_mean = post$sigma2_mean,
  #    lower = sigma2_lower,
  #    upper = sigma2_upper),
  #
  #  sigma = list(
  #    estimate_plugin = sigma_est_plugin,
  #    estimate_mean = sigma_est_mean,
  #    lower = sigma_lower,
  #    upper = sigma_upper))

  #print(npars)
  vec_lower <- vec_upper <- rep(NA, npars)

  idx_pis <- struct$pis$idx
  vec_lower[idx_pis] <- pi_lower
  vec_upper[idx_pis] <- pi_upper

  idx_mus <- struct$mus$idx
  vec_lower[idx_mus] <- as.vector(mu_lower)
  vec_upper[idx_mus] <- as.vector(mu_upper)

  # Be careful: lambda_as is defined K x D
  #   matrix(priors$a0, K, D, byrow = TRUE) +
  #   matrix(0.5 * colSums_rprobs, K, D, byrow = FALSE)
  # while struct uses component-order, so t(lambda_lower) is required


  idx_Lambdas <- unlist(lapply(struct$Lambdas, function(x) x$idx))
  vec_lower[idx_Lambdas] <- as.vector(t(lambda_lower))
  vec_upper[idx_Lambdas] <- as.vector(t(lambda_upper))
  #vec_lower[idx_Lambdas] <- as.vector(lambda_lower)
  #vec_upper[idx_Lambdas] <- as.vector(lambda_upper)

  res <- list(vec_lower = vec_lower, vec_upper = vec_upper)

  class(res) <- "mog_ng_cavi_ci"
  res
}

print.mog_ng_cavi_ci <- function(x, do_print = TRUE, digits = 4, parameter = c("lambda", "sigma2", "sigma"), ...)
{
  # TODO more informative version in 0.1.5,
  # here the output of confint.mog_normal_gamma_cavi was simplified
  # to stick to the requirements of the simulation

  #print(round(x$matrix, digits))
  cat("Lower bounds:\n")
  print(round(x$vec_lower, digits))
  cat("Upper bounds:\n")
  print(round(x$vec_upper, digits))

  #«# TODO create a separate function to build the matrix returned (invisibly) here
  #
  #parameter <- match.arg(parameter)
  #
  #level <- x$level
  #alpha <- 1 - level
  #
  #lower_name <- paste0(formatC(100 * alpha / 2, format = "fg"), " %")
  #upper_name <- paste0(formatC(100 * (1 - alpha / 2), format = "fg"), " %")
  #
  #make_ci_mat <- function(lower, upper) {
  #  out <- cbind(lower, upper)
  #  colnames(out) <- c(lower_name, upper_name)
  #  out
  #}
  #
  #K <- length(x$pis[, "estimate"])
  #D <- ncol(x$mus$estimate)
  #
  ## Mixture probabilities.
  #
  #pi_mat <- make_ci_mat(lower = x$pis[,"lower"], upper = x$pis[,"upper"])
  #rownames(pi_mat) <- paste0("pis[", seq_len(K), "]")
  #
  ## Component means.
  ##
  #
  #mu_rows <- expand.grid( k = seq_len(K), d = seq_len(D))
  #
  #mu_mat <- make_ci_mat(
  #  lower = x$mu$lower[cbind(mu_rows$k, mu_rows$d)],
  #  upper = x$mu$upper[cbind(mu_rows$k, mu_rows$d)])
  #
  #rownames(mu_mat) <- paste0("mus[", mu_rows$k, ",", mu_rows$d, "]")
  #
  ## Diagonal precision / variance / standard deviation parameters.
  #
  #if (parameter == "lambda") {
  #  par_lower <- x$Lambdas$lower
  #  par_upper <- x$Lambdas$upper
  #  par_prefix <- "Lambdas"
  #}
  #
  #if (parameter == "sigma2") {
  #  par_lower <- x$sigma2$lower
  #  par_upper <- x$sigma2$upper
  #  par_prefix <- "Sigmas2"
  #}
  #
  #if (parameter == "sigma") {
  #  par_lower <- x$sigma$lower
  #  par_upper <- x$sigma$upper
  #  par_prefix <- "Sigmas"
  #}
  #
  #par_rows <- expand.grid(
  #  k = seq_len(K),
  #  d = seq_len(D)
  #)
  #
  #par_mat <- make_ci_mat(
  #  lower = par_lower[cbind(par_rows$k, par_rows$d)],
  #  upper = par_upper[cbind(par_rows$k, par_rows$d)])
  #
  #rownames(par_mat) <- paste0(
  #  par_prefix, "[", par_rows$k, "][", par_rows$d, "]")
  #
  ## Final printed matrix.
  #
  #res <- rbind(pi_mat, mu_mat, par_mat)
  #
  #if (do_print)
  #  print(round(res, digits = digits), ...)
  #
  ##invisible(x)
  #invisible(res)
}


## -------------------------------------------------------------------
## Normal-Gamma family.
## -------------------------------------------------------------------

# Using update* naming convention would be neater, e.g.:
# update_auxiliary, update, update_responsibilities.
# However, using names not starting by the same substring can
# make list searches by name slightly easier/faster.

mog_normal_gamma <- list(
  update = mog_ng_upd,
  resp_update = mog_ng_upd_resp,
  init = mog_ng_init,
  elbo = mog_ng_elbo,
  posteriors = mog_ng_posteriors,
  name = "mog_normal_gamma"
)
