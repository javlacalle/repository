## Expectation-Maximization for a
## Vector Autoregressive multivariate mixture-of-Gaussians model.
## The cases or diagonal and full covariance matrices are considered.
##
## Model definition based on make_model().




# TODO
# A useful supplementary experiment would fit the inference network initially to the k-means responsibility matrix or to the exact initial EM responsibilities. Comparing that version with random initialization would show how much of the instability comes from the optimization landscape.

# NOTE learned from the simulation: Amortization is a parameterization, not automatically an inference guarantee.
#
#Its failures are informative rather than surprising: the network is heavily overparameterized, is randomly initialized, is not trained by an ELBO or posterior-matching criterion, and can steer the M-step toward singular mixture solutions. The simulation therefore supports the algorithm’s intended role as a conceptual experiment, but not as a competitive alternative to ordinary EM in this tractable setting.

# -------------------------------------------------------------------
# Expectation step: Update responsibilities.
# -------------------------------------------------------------------

em_log_scores <- function(model, include_constants = TRUE)
{
  # Debug:
  # if (determinant(pars$Lambdas[[k]], TRUE)$sign <= 0) stop("Matrix is not positive definite.")

  # NOTE documentation:
  # the scores are here the log-density values assigned to each observation at each specific component, 
  # aggregating the scores over all components gives the predicitive log-density (PDL) of each observation;
  # the PDL is numerically the same as the contribution of each observation to the likelihood 
  # (at least under when the observations are independent).
  # Conceptally they are not the same: PDL regards the observation as an outcome to be predicted,
  # while the likelihood contribution takes observations as observations to evaluation the likelihood function.

  if (model$is.diagonal)
  {
    log_scores_part <- function(E_k, Lambda_k)
      0.5 * sum(log(Lambda_k)) -
      # compare timings with: rowSums(sweep(E_k^2, 2L, Lambda_k, "*"))
      0.5 * colSums(Lambda_k * t(E_k^2))
  } else # model$is.fullcov
  {
    log_scores_part <- function(E_k, Lambda_k)
      0.5 * as.numeric(determinant(Lambda_k, TRUE)$modulus) -
      0.5 * rowSums((E_k %*% Lambda_k) * E_k)
  }

  X <- model$X
  pars <- model$pars
  N <- model$N
  K <- model$K

  Epart <- if (model$is.dynamic) X - model$Xlags %*% pars$B else X
  const <- if (include_constants) -0.5 * model$D * log(2 * pi) else 0.0

  log_scores <- matrix(NA, N, K)

  for (k in seq_len(K))
  {
    E_k <- sweep(Epart, 2L, pars$mus[k,], "-", FALSE) # Epart centered
    log_scores[,k] <- const +
      log(pars$pis[k]) + log_scores_part(E_k, pars$Lambdas[[k]])
  }

  log_scores
}

em_expectation <- function(model, debug = TRUE) # include_constants = FALSE
{
  # update responsibilities

  rprobs <- softmax_rows(em_log_scores(model, FALSE))

  if (debug && !isTRUE(all.equal(rowSums(rprobs), rep(1, model$N))))
    stop("Responsibilities do not sum to one.")

  rprobs
}


# -------------------------------------------------------------------
# Maximization step: Parameter updates.
# -------------------------------------------------------------------

em_maximization <- function(model)
{
  ridge_cov <- 1e-6

  X <- model$X
  Xlags <- model$Xlags
  K <- model$K
  N <- model$N
  D <- model$D

  rprobs <- model$rprobs
  colSums_rprobs <- colSums(rprobs)

  isDynamic <- model$is.dynamic
  isFullcov <- model$is.fullcov

  if (isDynamic) { E <- X - Xlags %*% model$pars$B } else E <- X

  fn_Sigmask_upd <-
    if (isFullcov) {
        function(k) crossprod((sqrt(rprobs[,k]) * E_centered_k)) / colSums_rprobs[k]
    } else # diagonal covariance matrix
        function(k) colSums(rprobs[,k] * E_centered_k^2) / colSums_rprobs[k]

  pis <- colSums_rprobs / N #sum(colSums_rprobs)
  #stopifnot(sum(colSums_rprobs) == N)

  mus <- (t(rprobs) %*% E) / colSums_rprobs

  # The diagonal entries of Lambdas_inv are the component- and
  # coordinate-specific weighted variances. Off-diagonal covariances
  # are constrained to zero.

  Lambdas_inv <- vector("list", K) # Sigmas

  for (k in seq_len(K))
  {
    E_centered_k <- sweep(E, 2L, mus[k,], "-", FALSE)
    Lambdas_inv[[k]] <- fn_Sigmask_upd(k)
  }

  # Lambdas <- lapply(Lambdas_inv, solve)
  # This may create numerical issues (singular matrix),
  # e.g. when a component has very few samples.
  # A ridge factor 1e-6 is added to avoid this.

  Lambdas <-
    if (isFullcov) {
      lapply(Lambdas_inv, function(x) solve(x + ridge_cov * diag(D)))
    } else
      lapply(Lambdas_inv, function(x) 1 / (x + ridge_cov))

  res <- list(pis = pis, mus = mus, Lambdas = Lambdas)

  if (isDynamic)
  {
    # ------------------------------------------------------------
    # Update B given current mus and Lambdas.
    #
    # Normal equation:
    #   sum_n sum_k tau_nk h_n h_n^T B Lambda_k =
    #   sum_n sum_k tau_nk h_n (x_n - mu_k)^T Lambda_k.
    #
    # Vectorized: M vec(B) = b.
    #
    # Using R column-major vectorization:
    #   vec(H B Lambda) = (Lambda^T \otimes H) vec(B).
    # ------------------------------------------------------------

    ridge_B <- 1e-8

    Dp <- D * model$p # ncol(model$Xlags)
    DpD <- Dp * D

    if (isFullcov) {
      Lambdas_mat <- Lambdas
    } else {
      # NOTE see if the calculations below can be simplified instead of converting here to matrix
      Lambdas_mat <- lapply(Lambdas, diag)
    }

    M_B <- matrix(0.0, DpD, DpD)
    b_B_mat <- matrix(0.0, Dp, D)

    for (k in seq_len(K))
    {
      Lambda_k <- Lambdas_mat[[k]]

      # Weighted cross-product sum_n tau_nk h_n h_n^T.
      H_weighted <- sqrt(rprobs[,k]) * Xlags
      S_hh_k <- crossprod(H_weighted)

      # Accumulate left-hand side matrix.
      M_B <- M_B + kronecker(Lambda_k, S_hh_k)

      # Accumulate right-hand side matrix:
      # sum_n tau_nk h_n (x_n - mu_k)^T Lambda_k.
      X_minus_mu_k <- sweep(X, 2L, mus[k,], "-", FALSE)

      b_B_mat <- b_B_mat + crossprod(Xlags, rprobs[,k] * X_minus_mu_k) %*% Lambda_k
    }

    M_B <- M_B + ridge_B * diag(DpD)
    B_vec <- solve(M_B, as.vector(b_B_mat))

    res$B <- matrix(B_vec, Dp, D)
  } # end (startsWith(model$family, "var_mog"))

  res
}


# -------------------------------------------------------------------
# Initialization of parameters in the EM procedure.
# -------------------------------------------------------------------

em_init_arcoefs <- function(model)
{
  # Initialize AR coefficients.
  # Used by em_init_resp() and pretrain_amortization_kmeans().

  if (model$is.dynamic)
  {
    # TODO allow pass these as arguments via '...'
    # also argument in em_var_mog_diagcov_upd()
    ridge_B <- 1e-8

    # Initialize B by ridge multivariate least squares.
    # Solve: B = (H^T H + ridge_B I)^(-1) H^T X,
    # where H is Xlags.

    X <- model$X
    Xlags <- model$Xlags
    Dp <- model$D * model$p # ncol(model$Xlags)

    XtX <- crossprod(Xlags)
    XtY <- crossprod(Xlags, X)
    B_init <- solve(XtX + ridge_B * diag(Dp), XtY)

    # VAR residuals: E = X - H B.
    # Xhat_init <- Xlags %*% B_init

    E <- X - Xlags %*% B_init

  } else { # !model$is.dynamic
    # stop("No AR coefficients to initialize; the input model is not dynamic.")
    # return the following instead of error message, in this way
    # if-else block is avoided in the functions that use this function.
    B_init <- NULL # required in the output
    E <- model$X
  }

  list(B = B_init, E = E)
}

em_init_resp <- function(model,
  method = c("kmeans", "random", "other"),
  kmeans_nstart = 10, kmeans_prob = 0.95,
  finit = NULL,  # used with method = "other" (currently no used, future flexibility)
  init_args = list(), # arguments passed to pretrain_amortization_kmeans() may be "added" here
  seed = NULL)
{
  # Initialize responsibilities
  #
  # NOTE documentation
  # For VAR-MoG models kmeans is applied on the residuals;
  # it requires initializing the matrix of coefficients,
  # which is then returned in B_init.
  #
  # NOTE documentation
  # method="other" used along with fn_init_resp and init_args;
  # currently not used, intended for future flexibility.
  #
  # init_args is used with method = "other" and with
  # "kmeans" when model$amortization$fn is "forward_network";
  # in the latter case it must define any of the arguments
  # to be passed to pretrain_amortization_kmeans().

  # The seed is used only when kmeans fails or with method "random".
  # The seed may therefore not be used.
  #
  # Alternatively to the block "verbose_seed" in mog_init(),
  # the simplest way to deal with the seed is to set it
  # once at the beginning (if it is not NULL) and reset to NULL
  # at the end of the function (regardless of whether it was used or not).

  # TODO see set default kmeans_nstart based on sample size

  if (!is.null(seed))
    set.seed(seed)

  # method <- match.arg(method) # comes already matched from em_init()

  N <- model$N
  K <- model$K

  am <- model$amortization

  # Initialize AR coefficients.
  #
  # Even if the model is not dynamic,
  # 'B_init' and 'E' need to be initialized.

  tmp <- em_init_arcoefs(model)
  B_init <- tmp$B
  E <- tmp$E

  if (model$is.amortized)
  {
    # Initialize etas: Array of dimension ncol(am$features) x K or a list
    # containing the coefficients of the inference function (amortization map).

    # TODO see whether the if-blocks below are better designed
    # splitting by argument 'method' instead of by is.null(am$fn)

    if (is.null(am$fn)) # "linear" amortization function by default
    {
      # By default apply a function linear in the parameters
      # (the features can be quadratic terms, ...).
      #
      # TODO linear in the parameters, so see do a linear regression:
      # etas <- rnorm(ncol(amortization$features) * K)
      # etas <- .lm.fit(x = amortization$features, y = X, tol = 1e-7)$coef
      # rprobs <- softmax_rows(as.vector(amortization$features %*% eta))

      if (method == "kmeans")
        stop(sprintf("Method '%s' is currently not implemented for linear amortization function. ", method),
             "Use method = 'random' or define the argument 'fn_init_resp' along with method='other'.")

      etas_nr <- NCOL(am$features)
      etas <- do.call("rnorm", c(list(n = etas_nr * (K - 1)), am$init_args))
      etas <- matrix(etas, etas_nr, K - 1)

    } else { # !is.null(am$fn)

      # Example: am$fn = "forward_network"; method = "random"

      #if (method == "other" && is.null(am$finit))
      #  stop("model$amortization$features is not NULL, but 'model$amortization$finit' is not defined ",
      #       "(use for example 'init_network_random').")

      if (method != "other")
      {
        if (am$fn != "forward_network")
          stop(sprintf("unexpected value for 'model$amortization$fn' = %s. ", am$fn),
            "Currently the methods 'kmeans' and 'random' support only 'forward_network'.")

        nn_width <- model$amortization$nn_width
        config_irandom <- list(input_dim = NCOL(am$features), width = nn_width, output_dim = 1)
      }

      if (method == "kmeans")
      {
        # NOTE documentation
        # kmeans initialization is not fully deterministic
        # in this setting (see comment below);
        # in fact, simulations show that it is advisable to try
        # different starting values.

        # preliminary random initialization:
        # these initial values are used to create the
        # slots in etas for the specified neural-network
        # (make_model() does not create those slots,
        # it just records the configuration of the network)
        #
        # these are also used as initial values for the
        # optim() procedure in pretrain_amortization_kmeans();
        # so they may have an effect on the result of avem_run()

        model$etas <- replicate(K - 1,
          do.call("init_network_random", config_irandom), simplify = FALSE)

        if (is.null(init_args$nstart)) init_args$nstart <- kmeans_nstart
        if (is.null(init_args$prob)) init_args$prob <- kmeans_prob
        init_args$full_output <- FALSE # overwrite if necessary, no problem

        etas <- do.call("pretrain_amortization_kmeans",
          c(list(model=model), init_args))

        # the following is done in em_init()
        # model$pars <- em_maximization(model)

      } else
      if (method == "random") {
        # Generate small random coefficients.

        #for (i in seq_len(K - 1)) do.call("init_network_random", am$args)
        etas <- replicate(K - 1,
          do.call("init_network_random", config_irandom), simplify = FALSE)

      } else { # method == "other"

        # NOTE Currently not used or tested; for future flexibility
        etas <- replicate(K - 1,
          do.call(fn_init_resp, init_args), simplify = FALSE)
      }

    } # end !is.null(am$fn)

    rprobs <- rprobs_from_amortization(etas, K, am, TRUE)

  } else # !model$is.amortized
  {
    # AR coefficients need to be initialized first to obtain
    # the residuals, upon which kmeans is applied

    rprobs <- NULL

    if (method == "kmeans")
    {
      km <- try(kmeans(E, centers = K, nstart = kmeans_nstart), silent = TRUE)

      if (inherits(km, "try-error")) {
        warning("K-means failed. Random initialization of responsibilities was applied.")
        # rprobs remains NULL
      } else {
        # K=1 is not expected.
        rprobs <- matrix((1 - kmeans_prob) / (K-1), N, K)
        rprobs[cbind(seq_len(N), km$cluster)] <- kmeans_prob
      }
    }

    if (method == "random" || is.null(rprobs))
    {
      #if (!is.null(seed))
      #  set.seed(seed)
      rprobs <- softmax_rows(matrix(rnorm(N * K), N, K))
    }

    # The parameters of the amortization function
    # need not be initialized in this case.
    etas <- NULL

  } # end if !model$is.amortized

  # complementary sanitation to that done in rprobs_from_amortization()
  if (any(colSums(rprobs) <= .Machine$double.eps))
  {
    warning("colSums(rprobs) <= .Machine$double.eps).",
      "\nSee note around this message in source code for a workaround to initialize 'pis' alternative to the update equation.")
  }

  if (!is.null(seed))
    set.seed(NULL)

  list(B = B_init, etas = etas, rprobs = rprobs)
}

em_init <- function(model,
  method = c("kmeans", "random", "other"), # for initialization of responsibilities
  kmeans_nstart = 10, kmeans_prob = 0.95, # for initialization of responsibilities
  fn_init_resp = NULL,  # used with method = "other" (currently no used, future flexibility)
  init_args = list(nn_width = 4), # passed to fn_init_resp or pretrain_amortization_kmeans()
  seed = NULL)
{
  # NOTE documentation:
  # 'seed' is used only if kmeans fails (if so, a warning is given); or with method "random".
  #
  # method = "other" used only to initialize responsibilities in a model with amortization.
  # fn_init_resp = NULL,  # used with method = "other" (currently no used, future flexibility)
  #
  # fn_init_resp and init_args are used to initialize the responsibilities in the amortized case;
  # the remaining parameters are initialized by means of the analytical expressions, so
  # no fn_init_resp or further arguments are needed to initialize them.
  # If method is kmeans or random there is no need to define fn_init_resp,
  # it is handled by default by em_init_resp.

  method <- match.arg(method)

  X <- model$X
  K <- model$K
  N <- model$N

  iresp <- em_init_resp(model, method, kmeans_nstart, kmeans_prob,
    fn_init_resp, init_args, seed)

  model$rprobs <- iresp$rprobs
  model$etas <- iresp$etas # NULL if is.null(model$amortization$features)
  if (!is.null(iresp$B))
  {
    # Do not do: if (model$is.dynamic):
    # If amortization is used or random initialization (instead of kmeans),
    # then it was not necessary to generate pars$B in em_init_resp();
    # if the model is dynamic it will be generated below by em_maximization().
    model$pars$B <- iresp$B
  }

  # Initialize the remaining parameters,
  # using the update equations obtained analytically in the maximization step.

  # update model parameters using analytical expressions from the maximization step
  #model$pars <- family$update(model)
  model$pars <- em_maximization(model) # ridge=1e-6

  # FIXME is this required? (in mog models, no var-mog, this is not done and seems reasonable)
  #
  # Recompute responsibilities using the full initialized model.
  # This makes rprobs consistent with pis, B, mus, and Lambdas.
  #
  #pars$rprobs <- em_var_mog_upd_resp(X = X, pars = pars,
  #  N = N, D = D, K = K,
  #  Xlags = Xlags,
  #  debug = TRUE)
  #
  # optional, one further update of the parameters from the latest responsibilities
  # new_pars <- em_var_mog_upd(...)
  # rprobs <- em_var_mog_upd_resp(newpars),

  model
}


# -------------------------------------------------------------------
# Complete-data observed log-likelihood and expected log-likelihood
# in the MoG-noise model.
# For monitoring convergence of EM, use the observed-data log-likelihood.
# -------------------------------------------------------------------

# em_mog_eloglik <- function(X, pars, debug = TRUE)
# {
#   # FIXME reuse code from em_mog_upd_resp() without duplicating it here.
#   # FIXME see if log_sum_exp required; also init eloglik <- 0 and  add all terms component by component.
#   K <- length(pars$pis)
#
#   # Include constants.
#
#   eloglik <- pars$rprobs * (-0.5 * ncol(X) * log(2*pi))
#
#   for (k in seq_len(K))
#     eloglik <- eloglik + sum(
#       pars$rprobs[,k] * (log(pars$pis[k]) + 0.5 * logdet_spd_mat(pars$Lambdas[[k]]) -
#       0.5 * quadratic_byrows(X, pars$mus[k,], pars$Lambdas[[k]])))
#
#   eloglik
# }

em_loglik <- function(model, include_constants = TRUE)
{
  log_scores <- em_log_scores(model, include_constants)

  sum(apply(log_scores, 1L, log_sum_exp))
}


# -------------------------------------------------------------------
# ELBO
# -------------------------------------------------------------------

em_elbo <- function(model, include_constants = TRUE)
{
  log_joint <- em_log_scores(model, include_constants)

  rprobs <- model$rprobs

  id <- rprobs > 0 # avoid log(0) in rprobs
  rprobs <- rprobs[id]

  sum(rprobs * (log_joint[id] - log(rprobs)))
}


# -------------------------------------------------------------------
# Coordinate-wise maximization of the
# expected complete-data log posterior in the VAR MoG-noise model.
# -------------------------------------------------------------------

em_run <- function(model,
  maxiter = 500, tol = 1e-8,
  verbose = FALSE, debug = FALSE)
  #do_sort = list(by_column = 1L, decreasing = FALSE)
{
  # Assumes proper definition of argument 'model',
  # generated by make_model() and initialized by em_init().

  K <- model$K
  N <- model$N
  D <- model$D

  # Storage.

  loglik_path <- rep(NA, maxiter + 1)
  loglik_path[1] <- em_loglik(model, FALSE)

  # Main loop.

  converged <- FALSE
  niter <- 0

  for (iter in seq_len(maxiter))
  {
    # Update responsibilities.

    # currently common function;
    # in the future more efficient version for diagonal covariance
    # may be defined as family$resp_update

    model$rprobs <- em_expectation(model, debug)

    # Update remaining parameters
    # using the analytical expressions from the maximization step

    #model$pars <- family$update(model)
    model$pars <- em_maximization(model) # ridge=1e-6

    # Check convergence and trace iterations.

    loglik_path[iter+1] <- em_loglik(model, FALSE)

    if (verbose)
    {
      cat(sprintf(
        "iter = %4d | logLik = %.6f | max mu = %.4f\n",
        iter, loglik_path[iter+1], max(abs(model$pars$mus))))
    }

    rel_change <- abs(loglik_path[iter+1] - loglik_path[iter]) / (1 + abs(loglik_path[iter]))

    niter <- iter

    #if (rel_change < tol && is.finite(rel_change))
    # This order avoids error if rel_change is NA.
    if (is.finite(rel_change) && rel_change < tol)
    {
      converged <- TRUE
      break
    }

  } # end main loop

  # Update responsibilities for the final updated parameters.
  # model$rprobs <- em_expectation(model, debug)

  sort.var_mog(model)

  res <- list(model = model, colSums_rprobs = colSums(model$rprobs),
    loglik_path = loglik_path[seq_len(niter + 1)],
    converged = converged, niter = niter)

  class(res) <- "em_fit"

  res
}


# -------------------------------------------------------------------
# Experimental.
# Hybrid EM in the VAR MoG-noise model:
#
#   o Numerical optimization of the parameters of an inference function
#     (linear function, neural-network) that parameterizes the
#     matrix of responsibilities.
#
#   o Analytical maximization step for the remaining parameters.
# -------------------------------------------------------------------

em_experimental_run <- function(model,
  optim_args = list(method = "BFGS", hessian = FALSE),
  do_sort = list(by_column = 1L, decreasing = FALSE),
  debug = FALSE)
{
  # Assumes proper elements in 'model' created by make_model()
  # and initialized by em_init().

  fn <- function(vec_etas)
  {
    model$etas <- vec_to_list_v2(vec_etas, etas_struct)

    model$rprobs <- rprobs_from_amortization(model$etas, K, model$amortization, FALSE)

    model$pars <- em_maximization(model)

    val <- em_loglik(model, FALSE)

    if (!is.finite(val)) return(.Machine$double.xmax)

    -val
  }

  K <- model$K

  etas <- model$etas
  etas_struct <- get_pars_structure_v2(etas)
  vec_etas <- if (is.list(etas)) unlist(etas, TRUE, FALSE) else as.vector(etas)

  opt <- do.call("optim", c(list(par = vec_etas, fn = fn), optim_args))

  # Final parameters and logLik evaluation.

  etas <- vec_to_list_v2(opt$par, etas_struct)
  rprobs <- rprobs_from_amortization(etas, K, model$amortization, debug)

  model$etas <- etas
  model$rprobs <- rprobs

  model$pars <- em_maximization(model)

  loglik <- em_loglik(model, FALSE)

  if (do_sort$by_column != 0)
  {

#FIXME adjust order when B AR coefs is present in var-mog models

#FIXME adjust etas of the amortization function (currently they are not used by other functions
# after the model is fitted, so at the moment it is not critical)

    # Assume proper arguments.
    #
    #if (!(do_sort$by_column %in% seq.int(nrow(mus))))
    #  stop("invalid value in 'do_sort$by_column'")
    #
    #if (is.null(do_sort$decreasing))
    #  do_sort$decreasing <- FALSE

    mus <- model$pars$mus
    ido <- order(mus[,do_sort$by_column], decreasing = do_sort$decreasing)

    if (!all(ido == seq.int(nrow(mus))))
    {
      model$pars$pis <- model$pars$pis[ido]
      model$pars$mus <- model$pars$mus[ido,,drop=FALSE]
      model$pars$Lambdas <- model$pars$Lambdas[ido]
      model$rprobs <- model$rprobs[,ido]
    }
  }

  res <- list(optim = opt, loglik = loglik, model = model)
  class(res) <- "em_experimental_fit"
  res
}

avem_run <- function(model,
  target = c("loglik", "elbo"),
  maxiter = 500, tol = 1e-8,
  optim_args = list(method = "BFGS", hessian = FALSE), # required only with amortize = TRUE
  verbose = FALSE, debug = FALSE)
  #do_sort = list(by_column = 1L, decreasing = FALSE)
{
  # Assumes proper definition of argument 'model',
  # generated by make_model() and initialized by em_init().

  # FIXME in future version remove argument 'target' and infer it from model definition;
  # for now keep it for experiments

  if (!model$is.classical) {
    cmod <- c('mog_full_cov',  'mog_diagonal_cov', 'var_mog_full_cov', 'var_mog_diagonal_cov')
    stop("avem_run() currently supports only model families:",
         paste(paste0("'", cmod, "'"), collapse = ","))
  }

  target <- match.arg(target)

  K <- model$K
  N <- model$N
  D <- model$D

  amortize <- model$is.amortized

  fn_target <- switch(target, "loglik" = em_loglik, "elbo" = em_elbo)

  if (amortize)
  {
    #etas <- model$etas
    etas_struct <- get_pars_structure_v2(model$etas)
    #vec_etas <- if (is.list(etas)) unlist(etas, TRUE, FALSE) else as.vector(etas)

    fn_a <- function(vec_etas)
    {
      model$etas <- vec_to_list_v2(vec_etas, etas_struct)

      model$rprobs <- rprobs_from_amortization(model$etas, K, model$amortization, FALSE)

      #model$pars <- em_maximization(model)

      #val <- em_loglik(model, FALSE)
      #val <- em_elbo(model, FALSE)
      val <- fn_target(model, FALSE)

      if (!is.finite(val)) return(.Machine$double.xmax)

      -val
    }

    # inference function for responsibilities
    fn_upd_rprobs <- function(model, debug)
    {
      # update responsibilities in the amortize scheme
      etas <- model$etas
      vec_etas <- if (is.list(etas)) unlist(etas, TRUE, FALSE) else as.vector(etas)

      opt <- do.call("optim", c(list(par = vec_etas, fn = fn_a), optim_args))

      etas <- vec_to_list_v2(opt$par, etas_struct)

      rprobs <- rprobs_from_amortization(etas, K, model$amortization, FALSE) # debug set to FALSE

      if (debug && !isTRUE(all.equal(rowSums(rprobs), rep(1, model$N))))
        stop("Responsibilities do not sum to one.")

      list(updates = list(etas = etas, rprobs = rprobs),
           optim = opt) # for further inspection of a fit
    }

  } else # !amortize
    fn_upd_rprobs <- function(model, debug)
      list(updates = list(rprobs = em_expectation(model, debug)),
           optim = list())

  # Storage.

  target_path <- rep(NA, maxiter + 1)
  target_path[1] <- fn_target(model, FALSE)
  optim_path <- vector("list", maxiter) # track optim results in amortized case

  # Main loop.

  converged <- FALSE
  niter <- 0

  for (iter in seq_len(maxiter))
  {
    # Expectation step. Update responsibilities.
    # update 'rprobs' and 'etas' in the amortized setting

    inference_upd <- fn_upd_rprobs(model, debug)
    model[names(inference_upd$updates)] <- inference_upd$updates

    optim_path[[iter]] <- inference_upd$optim

    # Maximization step. Update remaining parameters using
    # the analytical expressions from the maximization step.

    #model$pars <- family$update(model)
    model$pars <- em_maximization(model) # ridge=1e-6

    # Check convergence and trace iterations.

    target_path[iter+1] <- fn_target(model, FALSE)

    if (verbose)
    {
      # TODO see choose other than max(abs(model$pars$mus))
      cat(sprintf(
        "iter = %4d | %s = %.6f | max mu = %.4f\n",
        iter, target, target_path[iter+1], max(abs(model$pars$mus))))
    }

    rel_change <- abs(target_path[iter+1] - target_path[iter]) / (1 + abs(target_path[iter]))

    niter <- iter

    #if (rel_change < tol && is.finite(rel_change))
    # This order avoids error if rel_change is NA.
    if (is.finite(rel_change) && rel_change < tol)
    {
      converged <- TRUE
      break
    }

  } # end main loop

  # Update responsibilities for the final updated parameters.
  # model$rprobs <- em_expectation(model, debug)

  # Labels standarization, sort model parameters.
  # Options could be chosen as argument to this function. For now,
  # use default ordering (increasing by first column, see sort-model.R);
  # these options are not critical, omitting them
  # keeps the function interface simpler.

  model <- sort_var_mog(model)

  res <- list(target = target,
    target_path = target_path[seq_len(niter + 1)],
    optim_path = if (amortize) optim_path else NULL,
    model = model, colSums_rprobs = colSums(model$rprobs),
    converged = converged, niter = niter)

  class(res) <- "avem_fit"

  res
}
