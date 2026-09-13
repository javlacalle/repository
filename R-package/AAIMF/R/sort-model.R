
sort_var_mog <- function(x, by_column = 1L, decreasing = FALSE)
{
  # The sorting is now done outside avem_run() to simplify its code.
  # No major/relevant overhead is expected by doing so, for these reasons:
  # sort() is called only once after the iterations in avem_run().
  # Because of the copy-on-modify approach of R,
  # the unmodified elements in input 'x' are expected to keep their storage,
  # no copies or overhead is expected by returning here the full 'x'.
  # BTW Using '<<-' eg. x$pars$pis <- x$pars$pis[ido] does not work
  # because the object named 'x' is not defined in the scope of avem_run().

  # NOTE documentation:
  # In dynamic models, x$pars$B (AR coefficients) are currently shared by all mixture
  # components; so component-label sorting does not apply to them.

  if (!inherits(x, "var_mog"))
    stop("'model' is not of class 'var_mog'")

#FIXME adjust etas of the amortization function (currently they are not used by other functions
# after the model is fitted, so at the moment it is not critical)

  # Assume proper arguments.
  #
  #if (!(by_column %in% seq.int(nrow(mus))))
  #  stop("invalid value in 'by_column'")

  mus <- x$pars$mus
  ido <- order(mus[,by_column], decreasing = decreasing)

  if (all(ido == seq.int(nrow(mus))))
    return(x)

  # !all(ido == seq.int(nrow(mus)))
  x$pars$pis <- x$pars$pis[ido]
  x$pars$mus <- x$pars$mus[ido,,drop=FALSE] # use drop for safety
  x$pars$Lambdas <- x$pars$Lambdas[ido]
  x$rprobs <- x$rprobs[,ido,drop=FALSE] # as K>1, drop may not be necessary, keep it for safety

  if (x$K == 2 & is.list(x$etas))
  {
    # For K = 2, reordering of etas is simple:
    # the log-odds based on the amortization mapping are
    # log( q_etas (z=2 | data) / q_etas (z=1 | data) ).
    # Swapping the components labels implies a change of sign
    # log(q (z=2|x) / q (z=1|x)) = -log(q (z=1|x) / q (z=2|x)).
    # Therefore the suitable reordering seems to be a change of sign
    # in the coefficients.
    #
    # TODO do some checks to confirm, but seems correct.

    # Assumes proper definition of neural-network.
    # Assumes checks already done by make_model().

    etas_sorted <- x$etas[[1]]
    etas_sorted$W3 <- -etas_sorted$W3
    etas_sorted$b3 <- -etas_sorted$b3

    x$etas[[1]] <- etas_sorted
  }

  # TODO For K > 2 see the appropriate reordering; see also for the linear mapping.
  # Alternatively add argument do_sort = FALSE as in previous version in avem_run()

  msg <- c("\nThe parameters of the model and responsibilities where reordered ",
      "(to be explained in the documentation), ",
      "but the coefficients of the amortization mapping where not reorderd. ",
      "\nDo not use model$etas. Stick to use the responsibilities model$rprobs.")

  if (is.matrix(x$etas))
  {
    warning("Reordering of the coefficients of the amortization mapping is ",
      "currently implemented only for K = 2 and neural network amortization function.", msg)
  }

  if (x$K > 2)
  {
    warning("For K > 2 components, reordering of the coefficients of the amortization mapping ",
      "is not currently implemented.", msg)
  }

  x
}

# NOTE this is auxiliary for the simulations, rprobs are not ordered, just model pars;
#
# TODO define as model method; generalize to non-diagonal model, Normal-Wishart, ...

# TODO debug log-lik should not change after reordering
#
# loglik_before <- em_loglik(fit1$model)
# scores_before <- em_log_scores(fit1$model)
#
# tmp <- order_mog_diagcov_pars(
#   fit1$model$pars,
#   fit1$model$rprobs
# )
#
# fit1$model$pars <- tmp$opars
# fit1$model$rprobs <- tmp$orprobs
#
# loglik_after <- em_loglik(fit1$model)
# scores_after <- em_log_scores(fit1$model)
#
# stopifnot(isTRUE(all.equal(
#   loglik_before,
#   loglik_after,
#   tolerance = 1e-10
# )))
#
# stopifnot(isTRUE(all.equal(
#   scores_before[, tmp$order, drop = FALSE],
#   scores_after,
#   tolerance = 1e-10,
#   check.attributes = FALSE
# )))

# order_mog_diagcov_pars <- function(pars, rprobs, by_dim = 1L, decreasing = FALSE)
# {
#   # Assumes proper elements are defined in 'pars'.
#
#   #if (is.null(pars$mus))
#   #  stop("'pars' must contain 'mus'.")
#   #
#   #K <- nrow(pars$mus)
#   #D <- ncol(pars$mus)
#   #if (by_dim < 1L || by_dim > D)
#   #  stop("'by_dim' must be a valid column index of 'mu_ms'.")
#
#   ord <- order(pars$mus[,by_dim], decreasing = decreasing)
#
#   if (all(ord == seq.int(nrow(pars$mus)))) # no need to check length
#     return(list(opars = pars, orprobs = rprobs))
#
#   opars <- pars
#
#   if (!is.null(pars$pis))
#     opars$pis <- pars$pis[ord]
#
#   if (!is.null(pars$mus))
#     opars$mus <- pars$mus[ord,,drop=FALSE]
#
#   if (!is.null(pars$Lambdas))
#     opars$Lambdas <- pars$Lambdas[ord]
#     # fixed:
#     #for (i in seq_len(nrow(pars$mus)))
#     #  opars$Lambdas[[i]] <- pars$Lambdas[[i]][ord]
#
#   list(opars = opars, orprobs = rprobs[,ord])
# }

order_mog_ng_pars <- function(pars, rprobs, by_dim = 1L, decreasing = FALSE)
{
  # Assumes proper elements are defined in 'pars'.

  #if (is.null(pars$mu_ms))
  #  stop("'pars' must contain 'mu_ms'.")
  #
  #K <- nrow(pars$mu_ms)
  #D <- ncol(pars$mu_ms)
  #if (by_dim < 1L || by_dim > D)
  #  stop("'by_dim' must be a valid column index of 'mu_ms'.")

  ord <- order(pars$mu_ms[,by_dim], decreasing = decreasing)

  if (all(ord == seq.int(nrow(pars$mu_ms)))) # no need to check length
    return(list(opars = pars, orprobs = rprobs))

  opars <- pars

  if (!is.null(pars$pi_alphas))
    opars$pi_alphas <- pars$pi_alphas[ord]

  if (!is.null(pars$mu_ms))
    opars$mu_ms <- pars$mu_ms[ord,,drop=FALSE]

  if (!is.null(pars$mu_betas))
    opars$mu_betas <- pars$mu_betas[ord]

  if (!is.null(pars$lambda_as))
    opars$lambda_as <- pars$lambda_as[ord,,drop=FALSE]

  if (!is.null(pars$lambda_bs))
    opars$lambda_bs <- pars$lambda_bs[ord,,drop=FALSE]

  list(opars = opars, orprobs = rprobs[,ord])
}
