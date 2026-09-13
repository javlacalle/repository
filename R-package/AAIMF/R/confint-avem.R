# -------------------------------------------------------------------
# Confidence intervals for object of class 'avem_fit'.
# -------------------------------------------------------------------

# NOTE should pars$B be parameterized to fall in the region of stationarity, ... ???

pi_to_eta <- function(pis)
{
  # Baseline (pi[1] is the baseline) multinomial-logit parameterization.
  log(pis[-1L] / pis[1L])
}

eta_to_pi <- function(eta)
{
  exp_eta <- exp(eta)
  denom <- 1 + sum(exp_eta)
  c(1 / denom, exp_eta / denom)
}

chol_transform <- function(x)
{
  if (is.vector(x))
    return (log(x))

  # else is.matrix(x)
  R <- chol(x)
  diag(R) <- log(diag(R))

  # NOTE
  # No need to create the lower side, only the upper side will be used;
  # but be careful.
  #R[lower.tri(R, FALSE)] <- R[upper.tri(R, FALSE)]

  #R
  R[upper.tri(R, TRUE)]
}

chol_transform_undo <- function(x, shape = NULL, asmatrix = TRUE)
{
  # Assumes proper input 'shape' passed by caller.
  #
  # Converts a *vector* containing the *upper and diagonal elements*
  # of the transformed matrix, back to the original parametrization

  if (is.null(shape))
    return(exp(x))

  # else, build matrix
  # (actually shape[1] == shape[2], leave it this way to keep ot aligned with get_pars_structure())
  mat <- matrix(0, shape[1], shape[2])
  mat[upper.tri(mat, TRUE)] <- x
  # initializing mat with 0s instead of Nas avoids this explicitly
  #mat[lower.tri(mat, FALSE)] <- mat[upper.tri(mat, FALSE)]
  diag(mat) <- exp(diag(mat))
  mat <- crossprod(mat)

  if (asmatrix)
    return(mat)

  return(mat[upper.tri(mat, TRUE)])
}

#get_upper_mat <- function(mat) return(mat[upper.tri(mat, TRUE)])

confint.avem_fit <- function(object, level = 0.95, ...)
{
  if (!inherits(object, "avem_fit")) {
    stop("'object' must be an object of class 'avem_fit'.")
  }

  stopifnot(level >= 0 && level <= 1)

  n_warnings <- 0

  insert_na_at_pos <- function(x, idx_pi1)
    c(x[seq_len(idx_pi1 - 1L)], NA, x[seq(from = idx_pi1, to = length(x), by = 1L)])

  # TODO see keep only x as explicit argument in the definition
  fn <- function(x, model, struct, idx_pis, idx_Lambdas, shapes_Lambdas)
  {
    # recover the right size of 'x' to keep consistency with struct$idx
    # (recall that one of the probabilities was left out and recovered below)
    x <- insert_na_at_pos(x, idx_pis[1])

    #x[idx_Lambdas] <- exp(x[idx_Lambdas])
    # Actually the following would work for the models considered in this version
    # lapply(vs, function(z) chol_transform_undo(z, shape = struct$Lambdas[[1]]$shape))
    # Keep it anyway as below for generality, but I don't see a situation/model
    # where the dimensions of cov-matrices are not the same in all states.

    model$pars <- vec_to_list_v2(x, struct)

    #tmp <- plogis(x[idx_pis[-1]])
    #model$pars$pis <- c(1 - sum(tmp), tmp)
    exp_etas_pis <- exp(x[idx_pis[-1]])
    denom <- 1 + sum(exp_etas_pis)
    model$pars$pis <- c(1/denom, exp_etas_pis/denom)

    model$pars$Lambdas <- Map(function(idx, shp)
      chol_transform_undo(x[idx], shp), idx_Lambdas, shapes_Lambdas)

    -em_loglik(model)
  }

  relative_tol <- 1e-8

  model <- object$model
  pars <- model$pars
  #struct <- get_pars_structure(pars)
  struct <- get_pars_structure_v2(pars)
  idx_pis <- struct$pis$idx
  idx_Lambdas <- lapply(struct$Lambdas, function(x) x$idx)
  shapes_Lambdas <- lapply(struct$Lambdas, function(x) x$shape)

  z <- qnorm(1 - (1-level)/2)

  # transform parameters

  # Baseline multinomial-logit parametrisation.
  # pi[1] is the baseline.
  # The unconstrained vector has length K - 1 and stores
  # eta_k = log(pi_k / pi_1), k = 2, ..., K.

# FIXME if D>1 taking logs to off-diagonal elements in pars$Lambdas seems not appropriate
# also idx for pars$Lambdas should stick to diagonal and upper (or lower) elements, symmetric matrix

  tr_pars <- pars

  # this is valid for K=2; for K>2 may not ensure mixing probabilities add up to 1.
  #tr_pars$pis <- qlogis(pars$pis[-1])
  #
  # Baseline (pi[1] is the baseline) multinomial-logit parameterization.
  tr_pars$pis <- log(pars$pis[-1L] / pars$pis[1L])

  Lambdas <- pars$Lambdas
  if (is.vector(Lambdas[[1]])) {
    # NOTE current models use same shape for all Lambdas[[i]] for all i=1,...,K
    tr_pars$Lambdas <- lapply(Lambdas, function(x) log(x))
  } else {
    # these are *vectors* containing the diagonal and upper side elements
    # chol_transform_undo() convert it to a matrix in the original parametrization
    tr_pars$Lambdas <- lapply(Lambdas, chol_transform)
  }

  # a <- Lambdas[[1]]
  # vta <- chol_transform(a)
  # #b <- chol_transform_undo(at[upper.tri(at, TRUE)], struct$Lambdas[[1]]$shape)
  # b <- chol_transform_undo(vta, dim(a))
  # all.equal(a, b, check.attributes=FALSE)

  #mats <- lapply(struct$Lambdas, function(x) matrix(NA, x$shape[1], x$shape[2]))

  #v_tr_pars <- flatten_pars(tr_pars)
  v_tr_pars <- unlist(tr_pars, TRUE, FALSE)

  #H <- tryCatch(
  #  stats::optimHess(par = v_tr_pars, fn = fn, model = model, struct = struct),
  #  error = function(e) NULL)
  #if (is.null(H)) return(NULL)
  #
  # Keep the message (if any) returned optimHess().
  H <- stats::optimHess(par = v_tr_pars, fn = fn,
    model = model, struct = struct, idx_pis = idx_pis,
    idx_Lambdas = idx_Lambdas, shapes_Lambdas = shapes_Lambdas)

  # force symmetry
  H <- (H + t(H)) / 2

  evals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  if (any(evals <= 0)) {
    warning("Negative eigenvalues found. Observed information matrix H is not positive definite.", 
	"A local maximum of the log-likelihood may not have been reached.", 
	"Inverting H may not reflect the correct curvature. ")

    n_warnings <- 1
  }

  if (any(!is.finite(evals)) || min(evals) <= relative_tol * max(evals)) {
    # TODO see return here stop() instead of warning()
    warning("Non-finite or lower than threshold eigen-value(s) of the Hessian.")
    #"Observed information is not sufficiently positive definite."

    # eigenvalues are about linear combinations of parameters,
    # so the indices of troublesome eigen-values do not necessarily match indices of troublesome parameters;
    # the following preliminary approach is therefore not sensible/justified:
    #   max_eval <- max(evals)
    #   troublesome_rows <- !is.finite(evals) | (evals <= relative_tol * max_eval)
    #   Hid <- which(!troublesome_rows)
    #   solve(H[Hid,Hid])
    #
    # The following approach may be better to handle a troublesome observed-information matrix,
    # as it relies on eigen-directions rather than on parameters;
    # return anyway info on warning and discard these cases in the simulations
    max_eval <- max(evals)
    Hid <- is.finite(evals) & (evals > relative_tol * max_eval)
    evecs <- eigen(H, symmetric = TRUE, only.values = FALSE)$vectors
    Hinv <- evecs[,Hid,drop=FALSE] %*% diag(1/evals[Hid], nrow=sum(Hid)) %*% t(evecs[,Hid,drop=FALSE])

    n_warnings <- n_warnings + 1

  } else { # suitable observed-information matrix

    Hinv <- tryCatch(solve(H), error = function(e) NULL)
    #if (is.null(Hinv)) return(NULL)
    #Hinv <- solve(H)
  }

  se <- sqrt(pmax(diag(Hinv), 0))

  # compute the lower and upper bounds of the interval
  # for each set of parameters

  #y <- v_tr_pars
  #y_lower <- y - z * se
  #y_upper <- y + z * se

  y_lower <- v_tr_pars - z * se
  y_upper <- v_tr_pars + z * se

  y_lower <- insert_na_at_pos(y_lower, idx_pis[1])
  y_upper <- insert_na_at_pos(y_upper, idx_pis[1])

  # NOTE be careful,
  # because the inverse multinomial-logit map is nonlinear,
  # coordinatewise lower and upper endpoints may not preserve
  # the same order for all probabilities.

  eta_lower <- y_lower[idx_pis[-1]]
  eta_upper <- y_upper[idx_pis[-1]]
  pis_lower_raw <- eta_to_pi(eta_lower)
  pis_upper_raw <- eta_to_pi(eta_upper)
  y_lower[idx_pis] <- pmin(pis_lower_raw, pis_upper_raw)
  y_upper[idx_pis] <- pmax(pis_lower_raw, pis_upper_raw)

  # check these upon pis_x_raw, not y_lower[idx_pis], y_upper[idx_pis].
  # For K=2, lower and upper vectors often sum to one, because one probability is obtained as the complement of the other and the other.
  # But for K>2, the lower/upper bound need not sum to 1.
  stopifnot(all.equal(sum(pis_lower_raw), 1))
  stopifnot(all.equal(sum(pis_upper_raw), 1))

  #y_lower[idx_Lambdas] <- exp(y_lower[idx_Lambdas])
  #y_upper[idx_Lambdas] <- exp(y_upper[idx_Lambdas])
  Lambdas_lower <- Map(function(idx, shp)
    chol_transform_undo(y_lower[idx], shp, FALSE), idx_Lambdas, shapes_Lambdas)
  Lambdas_upper <- Map(function(idx, shp)
    chol_transform_undo(y_upper[idx], shp, FALSE), idx_Lambdas, shapes_Lambdas)

  idx_Lambdas <- unlist(idx_Lambdas, FALSE, FALSE)
  y_lower[idx_Lambdas] <- unlist(Lambdas_lower, FALSE, FALSE)
  y_upper[idx_Lambdas] <- unlist(Lambdas_upper, FALSE, FALSE)

  # output

  #list_lower <- vec_to_list_v2(y_lower, struct)
  #list_upper <- vec_to_list_v2(y_upper, struct)
  #vec_lower <- name_vec_pars(list_lower)
  #vec_upper <- name_vec_pars(list_upper)
  #
  #bounds_mat <- cbind(vec_lower, vec_upper)
  #colnames(bounds_mat) <- paste0(round(100*c((1-level)/2, (1+level)/2),3), " %")

  # NOTE despite 'bounds_mat' reports bounds for off-diagonal elements in
  # symmetric matrices (pars$Lambdas), the parameters passed to the
  # numerical Hessian include only the diagonal and upper side elements of the matrices.

  res <- list(
    #matrix = bounds_mat,
    vec_lower = y_lower, vec_upper = y_upper,
    #lower = list_lower, upper = list_upper,
    warnings = n_warnings)

  class(res) <- "em_ci"
  res
}

print.em_ci <- function(x, digits = 4, ...)
{
  #TODO see easier to interpret output (as former using matrix, but review names)

  if (!inherits(x, "em_ci")) {
    stop("'x' must be an object of class 'em_ci'.")
  }

  #print(round(x$matrix, digits))
  cat("Lower bounds:\n")
  print(round(x$vec_lower, digits))
  cat("Upper bounds:\n")
  print(round(x$vec_upper, digits))
}
