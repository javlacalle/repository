
# log-predictive densities computed in the simulations.

log_predictive_density <- function(model, X_test, B = 1000)
{
  if (!(model$family %in% c('mog_full_cov', 'mog_diagonal_cov', 'mog_normal_gamma', 'mog_normal_wishart')))
  	stop(sprintf("'log_predictive_density' is currently not implemented for model %s", model$family))

  # NOTE documentation
  # if model is 'mog_normal_gamma' or 'mog_normal_wishart', then 'test_pars'

  # NOTE
  # Modifying 'model_aux' as done below is in general risky:
  # model_aux$pars <- pars_b
  #
  # The slot 'pars' contains the priors, not the posteriors and
  # 'pars' will no longer be consistent the chosen 'priors'.
  #
  # As em_log_scores() uses only 'pars', 'model_aux' contains the right values.
  # This allows reusing em_log_scores() and simplifies the code.
  #
  # 'model_aux' is used only internally for the purposes mentioned before.
  # For now, this is the most simple approach to follow.

  model_aux <- model
  model_aux$X <- X_test

  if (model$family %in% c('mog_full_cov', 'mog_diagonal_cov'))
  {
    log_scores <- em_log_scores(model_aux, TRUE)
    return(apply(log_scores, 1L, log_sum_exp))
  }

  Ntest <- nrow(X_test)
  logdens_draws <- matrix(NA, Ntest, B)

  for (b in seq_len(B))
  {
    pars_b <- var_mog_draw_pars(model)

    # old version
    #logdens_draws[,b] <- dmog_diag_logdens(X = Xtest,
    #  pis = pars_b$pis, mus = pars_b$mus, Lambdas = pars_b$Lambdas)

    model_aux$pars <- pars_b
    logdens_draws[,b] <- apply(em_log_scores(model_aux, TRUE), 1L, log_sum_exp)
  }

  # For each test observation:
  # log E_q[p(x_new | theta)]
  # approx log { B^{-1} sum_b p(x_new | theta_b) }.

  vapply(seq_len(Ntest),
    function(n) log_sum_exp(logdens_draws[n,]) - log(B), numeric(1))
}

var_mog_draw_pars <- function(model)
{
  # TODO Draw parameters from the mog_normal_wishart() (see also var_mog models)

  if (!(model$family %in% c('mog_normal_gamma')))
    stop(sprintf("'var_mog_draw_pars' is currently not implemented for model %s", model$family))

  pars <- model$pars
  K <- model$K
  D <- model$D

  pis <- as.numeric(gtools::rdirichlet(1L, pars$pi_alphas))

  mus <- matrix(NA, nrow = K, ncol = D)
  lambda_mat <- matrix(NA, nrow = K, ncol = D)

  for (k in seq_len(K))
  {
    for (d in seq_len(D))
    {
      lambda_mat[k, d] <- rgamma(1, shape = pars$lambda_as[k, d], rate = pars$lambda_bs[k, d])

      mus[k, d] <- rnorm(1, mean = pars$mu_ms[k, d],
        sd = sqrt(1 / (pars$mu_betas[k] * lambda_mat[k, d])))
    }
  }

  Lambdas <- lapply(seq_len(K), function(k) lambda_mat[k, ])

  list(pis = pis, mus = mus, Lambdas = Lambdas)
}
