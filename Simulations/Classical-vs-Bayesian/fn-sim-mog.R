
# auxiliary functions

# is_inside <- function(v, x) x >= v[1] && x <= v[2]
#
# #FIXME use  ci$bounds_mat
# are_inside <- function(ci, dgp)
#   c(is_inside(ci$bounds_mat["pis[1]",], dgp$pis[1]),
#     is_inside(ci$bounds_mat["pis[2]",], dgp$pis[2]),
#     is_inside(ci$bounds_mat["mus[1,1]",], dgp$mus[1,1]),
#     is_inside(ci$bounds_mat["mus[2,1]",], dgp$mus[2,1]),
#     is_inside(ci$bounds_mat["mus[1,2]",], dgp$mus[1,2]),
#     is_inside(ci$bounds_mat["mus[2,2]",], dgp$mus[2,2]),
#     is_inside(ci$bounds_mat["Lambdas[1][1]",], dgp$Lambdas[[1]][1]),
#     is_inside(ci$bounds_mat["Lambdas[1][2]",], dgp$Lambdas[[1]][2]),
#     is_inside(ci$bounds_mat["Lambdas[2][1]",], dgp$Lambdas[[2]][1]),
#     is_inside(ci$bounds_mat["Lambdas[2][2]",], dgp$Lambdas[[2]][2]))

ci_get_info <- function(ci, struct_mod1, dgp_vec_pars, dgp_struct) #debug = TRUE)
{
  # This version is safer as it avoids issues with possible different orderings
  # when packing and unpacking pars in the fitted model and in DGP pars.

  # debug: omit in simulations
  #if (debug)
  #{
  #  struct_lower <- get_pars_structure_v2(ci$lower)
  #  struct_upper <- get_pars_structure_v2(ci$upper)
  #  stopifnot(identical(struct_lower, struct_mod1))
  #  stopifnot(identical(struct_upper, struct_mod1))
  #}

  #if (is.list(ci$lower))
  #{
  #  ci_vec_lower <- unlist(ci$lower)
  #  ci_vec_upper <- unlist(ci$upper)
  #} else {
  #  # 'ci' based on 'mog_cavi' returns lower and upper already as vectors
  #  ci_vec_lower <- ci$lower
  #  ci_vec_upper <- ci$upper
  #}

  ci_vec_lower <- ci$vec_lower
  ci_vec_upper <- ci$vec_upper

  ci_coverage <- dgp_vec_pars * NA

  for (ref in setdiff(names(struct_mod1), "Lambdas"))
  {
    id1 <- dgp_struct[[ref]]$idx
    id2 <- struct_mod1[[ref]]$idx

    bmat <- cbind(ci_vec_lower[id2], ci_vec_upper[id2], dgp_vec_pars[id1])
    ci_coverage[id1] <- apply(bmat, 1, function(x) x[3] >= x[1] && x[3] <= x[2])
  }

  for (k in seq_along(dgp_struct[["Lambdas"]]))
  {
    id1 <- dgp_struct[["Lambdas"]][[k]]$idx
    id2 <- struct_mod1[["Lambdas"]][[k]]$idx

    bmat <- cbind(ci_vec_lower[id2], ci_vec_upper[id2], dgp_vec_pars[id1])
    ci_coverage[id1] <- apply(bmat, 1, function(x) x[3] >= x[1] && x[3] <= x[2])
  }

  list(coverage = as.logical(ci_coverage), width = ci_vec_upper - ci_vec_lower)
}

#diff_lists <- function(a, b)
#{
#  #stopifnot(is.list(a), is.list(b))
#  #stopifnot(identical(names(a), names(b)))
#
#  diff_elem <- function(x, y)
#  {
#    if (is.list(x)) {
#      # same structure assumed for y
#      out <- vector("list", length(x))
#      #names(out) <- names(x)
#      #for (nm in names(x)) out[[nm]] <- diff_elem(x[[nm]], y[[nm]])
#      for (i in seq_along(x)) out[[i]] <- diff_elem(x[[i]], y[[i]])
#      out
#    } else {
#      # vectors/matrices/arrays: subtract directly
#      x - y
#    }
#  }#
#
#  b_names <- names(b)
#  out <- vector("list", length(a))
#  names(out) <- names(a)
#  for (nm in names(a))
#    if (nm %in% b_names)
#      out[[nm]] <- diff_elem(a[[nm]], b[[nm]])
#  out
#}

mog_ng_posterior_probs <- function(model, B = 5000)
{
  #if (is.null(model$pars))
  #  stop("'model' must contain 'pars'.")

  pars <- model$pars
  K <- model$K
  D <- model$D

  if (K != 2L)
    stop("This helper is currently written for K = 2.")

  # ------------------------------------------------------------
  # Posterior probability for pi_2 > pi_1.
  #
  # q(pi) = Dirichlet(alpha_1, alpha_2).
  # For K = 2, pi_2 marginally follows:
  #
  # pi_2 ~ Beta(alpha_2, alpha_1).
  #
  #
  # Pr(pi_2 > pi_1) = Pr(pi_2 > 0.5).
  # ------------------------------------------------------------

  alpha1 <- pars$pi_alphas[1L]
  alpha2 <- pars$pi_alphas[2L]

  prob_pi2_gt_pi1 <- 1 - pbeta(0.5, shape1 = alpha2, shape2 = alpha1)

  # ------------------------------------------------------------
  # Posterior probabilities for mean comparisons.
  #
  # For the Normal-Gamma variational posterior:
  #
  # lambda_kd ~ Gamma(a_kd, b_kd)
  # mu_kd | lambda_kd ~ N(m_kd, (beta_k lambda_kd)^(-1)).
  #
  # The marginal distribution of mu_kd is Student-t.
  # For simplicity the probability of mu_21 > mu_11 is
  # approximated here by simulation from mu_kd | lambda_kd .
  # ------------------------------------------------------------

  mu_draws <- array(NA, dim = c(B, K, D))

  for (k in seq_len(K))
  {
    for (d in seq_len(D))
    {
      lambda_draws <- rgamma(B, shape = pars$lambda_as[k, d], rate = pars$lambda_bs[k, d])
      mu_draws[, k, d] <- rnorm(B,
        mean = pars$mu_ms[k, d], sd = sqrt(1 / (pars$mu_betas[k] * lambda_draws)))
    }
  }

  prob_mu21_gt_mu11 <- mean(mu_draws[,2L,1L] > mu_draws[,1L,1L])

  prob_mu22_gt_mu12 <- mean(mu_draws[,2L,2L] > mu_draws[,1L,2L])

  c(VB_Pr_pi2_gt_pi1 = prob_pi2_gt_pi1,
    VB_Pr_mu21_gt_mu11 = prob_mu21_gt_mu11,
    VB_Pr_mu22_gt_mu12 = prob_mu22_gt_mu12)
}

# older version (an updated version is currently included in the AAIMF package):

#logsumexp <- function(x)
#{
#  m <- max(x)
#  m + log(sum(exp(x - m)))
#}

#dmog_diag_logdens <- function(X, pis, mus, Lambdas)
#{
#  # This function could be replaced by em_log_scores(), 
#  # available in the latest version of the AAIMF package
#  # em_log_scores() is more general interface that handles different models.
#  # For now keep this previous version to compute the predictive densities 
#  # in the simulations.
#
#  N <- nrow(X)
#  K <- length(pis)
#  D <- ncol(X)
#
#  hDlog2pi <- 0.5 * D * log(2 * pi)
#
#  out <- numeric(N)#
#
#  for (n in seq_len(N))
#  {
#    log_terms <- numeric(K)
#
#    for (k in seq_len(K))
#    {
#      lambda_k <- Lambdas[[k]]#
#
#      x_centered <- X[n,] - mus[k,]
#
#      log_terms[k] <- log(pis[k]) + 0.5 * sum(log(lambda_k)) -
#        hDlog2pi - 0.5 * sum(lambda_k * x_centered^2)
#    }
#
#    out[n] <- logsumexp(log_terms)
#  }
#
#  out
#}

#mog_ng_draw_pars <- function(model)
#{
#  pars <- model$pars
#  K <- model$K
#  D <- model$D
#
#  pis <- as.numeric(gtools::rdirichlet(1L, pars$pi_alphas))
#
#  mus <- matrix(NA, nrow = K, ncol = D)
#  lambda_mat <- matrix(NA, nrow = K, ncol = D)
#
#  for (k in seq_len(K))
#  {
#    for (d in seq_len(D))
#    {
#      lambda_mat[k, d] <- rgamma(1, shape = pars$lambda_as[k, d], rate = pars$lambda_bs[k, d])
#
#      mus[k, d] <- rnorm(1, mean = pars$mu_ms[k, d],
#        sd = sqrt(1 / (pars$mu_betas[k] * lambda_mat[k, d])))
#    }
#  }
#
#  Lambdas <- lapply(seq_len(K), function(k) lambda_mat[k, ])
#
#  list(pis = pis, mus = mus, Lambdas = Lambdas)
#}

#mog_ng_posterior_predictive_logdens <- function(Xtest, pis, mus, Lambdas, B = 1000)
#{
#  # TODO Adjust this function to use the more general interface em_log_scores()
#  # instead of dmog_diag_logdens().
#  # It will require modifying the slot 'pars' in the model;
#  # not so straightforward because they should be kept consistent with the chosen priors;
#  # but internally, for the purposes here (call em_logscores()) it could be
#  # done just by model$pars$mus <- mus, ...
#
#  Ntest <- nrow(Xtest)
#  logdens_draws <- matrix(NA, Ntest, B)
#
#  for (b in seq_len(B))
#  {
#    pars_b <- mog_ng_draw_pars(model)
#
#    # old version
#    logdens_draws[,b] <- dmog_diag_logdens(X = Xtest,
#      pis = pars_b$pis, mus = pars_b$mus, Lambdas = pars_b$Lambdas)
#  }
#
#  # For each test observation:
#  #
#  # log E_q[p(x_new | theta)]
#  # approx log { B^{-1} sum_b p(x_new | theta_b) }.
#
#  vapply(seq_len(Ntest),
#    function(n) logsumexp(logdens_draws[n,]) - log(B), numeric(1))
#}
