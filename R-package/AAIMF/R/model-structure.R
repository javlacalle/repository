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

get_pars_structure <- function(pars)
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

  build <- function(x)
  {
    if (is.list(x))
      return(lapply(x, build))

    d <- dim(x)
    len <- if (is.null(d)) length(x) else prod(d)
    idx <- seq.int(id, length.out = len)
    id <<- id + len # scope within get_pars_structure().

    list(shape = d, length = len, idx = idx)
  }

  build(pars)
}

vec_to_list <- function(vec_pars, struct, transform = list(NULL, NULL))
{
  # Convert the flattened vector of parameters to a list as defined in 'struct'.
  # Recursive calls made in turn by lapply().

  if (!is.null(transform[[1]]))
  {
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

  fill_as_struct <- function(s)
  {
    if (all(c("shape", "idx") %in% names(s)))
    {
      vals <- vec_pars[s$idx]
      dim(vals) <- s$shape
      return(vals)
    } #else: Skip, the element is not the innermost list
      # containing 'shape' and 'idx', e.g. 'Lambda_Ws'.

    lapply(s, fill_as_struct)
  }

  #FIXME See set gammas to NULL instead of list()

  fill_as_struct(struct)
}

#flatten_pars <- function(pars) unlist(pars, use.names = FALSE)

name_vec_pars <- function(x)
{
  make_names <- function(vec, parent)
  {
    n <- length(vec)
    if (is.null(dim(vec))) {
      # vector
      idx <- seq_len(n)
      paste0(parent, "[", idx, "]")
    } else {
      # array (here: matrices), use arrayInd to get [i,j]
      d <- dim(vec)
      idx_mat <- arrayInd(seq_len(n), .dim = d)
      # build [i,j] names
      apply(idx_mat, 1, function(ij) {
        paste0(parent, "[", ij[1], ",", ij[2], "]")
      })
    }
  }

  nm <- character(0)

  for (parent_name in names(x)) {
    obj <- x[[parent_name]]

    if (is.list(obj) && !is.null(names(obj))) {
      # not needed for your example, but kept generic
      for (k in seq_along(obj)) {
        if (is.null(dim(obj[[k]]))) {
          nm <- c(nm, make_names(obj[[k]], paste0(parent_name, "[", k, "]")))
        } else {
          nm <- c(nm, make_names(obj[[k]], paste0(parent_name, "[", k, "]")))
        }
      }
    } else if (is.list(obj)) {
      # list-of-vectors (your Lambdas)
      for (k in seq_along(obj)) {
        nm <- c(nm, make_names(obj[[k]], paste0(parent_name, "[", k, "]")))
      }
    } else {
      # vector or matrix/array (your pis, mus)
      nm <- c(nm, make_names(obj, parent_name))
    }
  }

  # assign to unlist
  u <- unlist(x, recursive = TRUE, use.names = FALSE)
  names(u) <- nm
  u
}

# x <- list(
#   pis = c(0.43, 0.57),
#   mus = structure(c(-1.88, 2.02, -0.93, 0.98), dim = c(2L, 2L)),
#   Lambdas = list(c(3.02, 1.45), c(0.98, 1.58))
# )
#
# u <- name_unlist(x)
# u
# names(u)
