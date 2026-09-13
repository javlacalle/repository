dgp_sanity_checks <- function(dgp, model_family)
{
  # Create auxiliary elements and do sanity checks.

  # Further auxiliary elements.
  Lambdas <- switch(model_family,
    "mog_diagonal_cov" = lapply(dgp$Sigmas, function(x) 1/diag(x)),
    "mog_full_cov" = lapply(dgp$Sigmas, function(x) solve(x)))

  dgp$Lambdas <- Lambdas
  dgp_pars_nms <- setdiff(names(dgp), c("N", "Sigmas"))
  dgp_vec_pars <- unlist(dgp[dgp_pars_nms])

  # Sanity check.
  if (!all(order(dgp$mus[,1], decreasing=FALSE) == seq.int(nrow(dgp$mus)))) {
    stop("Error: 'mus[,1]' must be defined in increasing order; for compatibility with pars order below.")
  } else
    cat(paste("OK DGP definition.\n"))

  # Sanity check:
  # Relevant for point_bias within one_iteration() and for the summary of results.
  #
  # Ensure unlist() flattens the list of parameters
  # into a vector containing the parameters in the same order;
  # if an error is get, it is because the elements defined in the list 'dgp'
  # at the beginning of this file are not defined in the same order as those
  # in fit$model$pars (by default "pis", then "mus", then "Lambdas").

  dgp_data <- mog_sim(N = 50, pis = dgp$pis, mus = dgp$mus, Sigmas = dgp$Sigmas)
  model0 <- em_init(make_model(model_family, dgp_data$X, dgp_data$K))
  isok <- identical(names(dgp_vec_pars), names(unlist(model0$pars)))
  msg <- "The order of the parameters in the DGP specification and in the model match each other.\n"
  if (isok) {
    cat(paste("OK", msg))
  } else
    stop(paste("Error:", sub("match", "do not match", msg)))

  #rm(dgp_data, model0, isok, msg)
  list(Lambdas = Lambdas, 
    dgp_pars_nms = dgp_pars_nms, dgp_vec_pars = dgp_vec_pars,
    dgp_struct = get_pars_structure_v2(dgp[dgp_pars_nms]))
}
