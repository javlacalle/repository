# -------------------------------------------------------------------
# Model structure.
# Functions to handle the management of model parameters depending
# on the needs (eg. stats::optim() requires a vector of parameters).
#
# The parameters of the model are defined by type in a list
# (see mog_nw_init(), mog_ng_init()).
#
# The list can be flattened to a raw vector, as required by stats::optim().
# Flattened vectors, matrices and nested lists (neural-network parameters),
# are converted back to an easier to read list of parameters according to
# the structure obtained by get_pars_structure().
# -------------------------------------------------------------------

# NOTE
# softplus: \( \tau = \log\big(1+e^{\eta}\big) \)
# Consider softplus transformation to impose the non-negative constraint on a variance (precision) parameter:
# it seems it allows the parameter o be able to get very close to \(0\) more stably

# NOTE documentation (validation example):
# struct <- get_pars_structure(pars)
# reconstructed <- vec_to_list(flatten_pars(pars), struct)
# isTRUE(all.equal(flatten_pars(reconstructed), flatten_pars(pars)))

get_pars_structure_v2 <- function(pars)
{
  # Recursive calls.
  # By applying lapply() to all the elements in 'pars' and by
  # recursive calls to build(), the recursions get access to all
  # sublists with each element in 'pars'.
  # This approach is cleaner than previous versions and can
  # handle the sublists 'Lambda_Ws'.
  #
  # Returns a nested list with the same elements as 'pars',
  # but instead of containing the values in 'pars', the shape,
  # length and indices of each element are given as output.
  # The indices refer to the location in a flattened version
  # of 'pars' (unlist(), flatten_pars() below).

  id <- 1L

  build_general <- function(x)
  {
    if (is.list(x))
      return(lapply(x, build_general))

    d <- dim(x)
    len <- if (is.null(d)) length(x) else prod(d)
    idx <- seq.int(id, length.out = len)
    id <<- id + len # scope within get_pars_structure().

    list(shape = d, length = len, idx = idx)
  }

  if (is.null(names(pars))) # etas parameters of amortization neural-network
    return(build_general(pars))

  build_symmetric_mat <- function(x)
  {
    if (is.list(x))
      return(lapply(x, build_symmetric_mat))

    d <- dim(x)
    stopifnot(length(d) == 2 && isTRUE(all(d[1] == d[2])))
    len <- if (is.null(d)) length(x) else d[1] * (d[1] + 1) / 2
    idx <- seq.int(id, length.out = len)
    id <<- id + len # scope within get_pars_structure().

    list(shape = d, length = len, idx = idx)
  }

  res <- vector("list", length(pars))
  names(res) <- names(pars)

  # NOTE documentation:
  # There are also matrices that are not vcov matrices,
  # so names must be used to decide how to deal with matrices;
  # be careful when extending the code an giving names to new parameters.
  mask <- names(pars) %in% c("Lambdas", "Sigmas")

  # adjust 'mask' for model with diagonal precision/vcov matrix instead of full matrix
  if (sum(mask) > 1)
    stop("TODO check case 'sum(mask) > 1'.")
  mask_mats <- unique(unlist(lapply(pars[mask][[1]], function(x) is.matrix(x))))
  if (length(mask_mats) != 1)
    stop("TOCHECK unexpected 'length(mask_mats) != 1'")
  mask[mask] <- mask_mats

  for (i in seq_along(pars))
    res[[i]] <- if (mask[i]) build_symmetric_mat(pars[[i]]) else build_general(pars[[i]])

  res
}

vec_to_list_v2 <- function(vec_pars, struct, transform = list(NULL, NULL))
{
  # Convert the flattened vector of parameters to a list as defined in 'struct'.
  # Recursive calls made in turn by lapply().

  if (!is.null(transform[[1]]))
  {
    # FIXME see do transformations in confint() here

    #names <- transform[[1]]
    # NOTE Do not do: fn <- match.fun(transform[[2]][i]), because of
    #Error in transform[[2]][i] : object of type 'builtin' is not subsettable.

    fns <- transform[[2]]
    if (is.null(fns))
      stop("transform[[2]] must be provided")

    # allow a single function to apply to all names
    if (!is.list(fns)) fns <- rep(list(fns), length(transform[[1]]))
      stopifnot(length(transform[[1]]) == length(fns))

    for (i in seq_along(transform[[1]]))
    {
      idx <- struct[[ transform[[1]][[i]] ]]$idx
      fn <- match.fun(fns[[i]])
      vec_pars[idx] <- fn(vec_pars[idx])
    }
  }

  fill_general <- function(s)
  {
    if (all(c("shape", "idx") %in% names(s)))
    {
      vals <- vec_pars[s$idx]
      dim(vals) <- s$shape
      return(vals)
    } #else: Skip, the element is not the innermost list
      # containing 'shape' and 'idx', e.g. 'Lambda_Ws'.

    lapply(s, fill_general)
  }

  if (is.null(names(struct))) # etas parameters of amortization neural-network
    return(fill_general(struct))

  fill_symmetric_mat <- function(s)
  {
    if (all(c("shape", "idx") %in% names(s)))
    {
        mat <- matrix(0, s$shape[1], s$shape[2])
        mat[upper.tri(mat, TRUE)] <- vec_pars[s$idx]
        mat[lower.tri(mat, FALSE)] <- mat[upper.tri(mat, FALSE)]
        return(mat)
    }
    lapply(s, fill_symmetric_mat)
  }

  res <- vector("list", length(struct))
  names(res) <- names(struct)

  # NOTE documentation:
  # There are also matrices that are not vcov matrices,
  # so names must be used to decide how to deal with matrices;
  # be careful when extending the code an giving names to new parameters.
  #
  #mask <- names(res) %in% c("Lambdas", "Sigmas")
  mask <- names(res) %in% "Lambdas"

  # adjust 'mask' for model with diagonal precision/vcov matrix instead of full matrix
  if (sum(mask) > 1)
    stop("TODO check case 'sum(mask) > 1'.")
  #mask_mats <- unique(unlist(lapply(struct[mask], function(x) !is.null(x$shape))))
  mask_mats <- unique(unlist(lapply(struct$Lambdas, function(x) !is.null(x$shape))))
  if (length(mask_mats) != 1)
    stop("TOCHECK unexpected 'length(mask_mats) != 1'")
  mask[mask] <- mask_mats

  for (i in seq_along(struct))
    res[[i]] <- if (mask[i]) fill_symmetric_mat(struct[[i]]) else fill_general(struct[[i]])

  res
}
