## -------------------------------------------------------------------
## Compute responsibilities inferred by the amortization function.
## -------------------------------------------------------------------

rprobs_from_amortization <- function(etas, K,
  amortization = list(features = NULL, fn = NULL, args_fn = list()),
  debug = TRUE)
{
  # NOTE The arguments are collected in a list for convenience;
  # this helps organizing the arguments of other functions (mog_init, mog_elbo).

  # NOTE documentation:
  # Currently 'amortization' may contain just the element 'features';
  # depending on the class of 'etas' matrix/list, by default a
  # linear/neural-network amortization is used.

  # NOTE documentation:
  # If etas is a list of parameters (neural-network amortization function)
  # it should contain K-1 elements.
  # In this case the argument K is ignored.
  #
  # If etas is a matrix (linear amortization function) and
  # when K > 2, then it arises the question whether it
  # is a reasonable approach completing with multiple columns of
  # zeros before applying softmax_rows (see comment below).

  # Compared to calling several times to amortization$arg, extracting the
  # arguments once at the beginning may reduce overhead (probably negligible).
  # It also reduces the length of the syntax.

  features <- amortization$features
  fn <- amortization$fn
  args <- amortization$args_fn

  if (is.matrix(etas))
  {
    if (is.null(fn))
    {
      # By default use a linear amortization function.
      # Linear in the parameters, a quadratic polynomial is allowed,
      # eg. features = cbind(1, x, x^2).

      fn <- function(X, etas) X %*% etas
    }

    # previous version (assumed K=2)
    # Metas <- cbind(0, etas)
    # rprobs <- softmax_rows(do.call(fn,
    #     c(list(X = features, etas = Metas), args)))

    # evaluate the amortization function

    pseudo_rprobs <- do.call(fn, c(list(X = features, etas = etas), args))

    # apply softmax_rows() to convert 'pseudo_rprobs' to probabilities
    #
    # 'pseudo_rprobs' may need to be completed with column(s) of zeros,
    # so that 'rprobs' contains K columns.
    #
    # Example: in the linear case 'pseudo_rprobs' is a column vector,
    # applying softmax_rows() to it would transform it a column of 1s,
    # if a column of zeros is appended, K=2 probabilities summing to 1 are
    # returned by softmax_rows() for each observation (N x K).
    #
    # If fn() returns K columns, then zero_mat has 0 columns,
    # cbind(zero_mat, pseudo_rprobs) is allowed (no concatenation).
    #
    # FIXME If K > 2, preappending 2 columns of zeros returns also probabilities,
    # but is this correct in this case, any issue when multiple columns are zeros ???

    zero_mat <- matrix(0, nrow(features), K - NCOL(pseudo_rprobs))
    rprobs <- softmax_rows(cbind(zero_mat, pseudo_rprobs))

  } else
  if (is.list(etas))
  {
    # Neural-network amortization function.
    #
    # Functions other than a neural-network that require
    # defining arguments as a list are also allowed.

    if (is.null(fn))
    {
      # By default a neural-network amortization function is used.
      # "forward_network" is defined in the source file neural-network.R of the package.

      fn <- "forward_network"
      #args <- NULL # warning if amortization$args_fn is not NULL (unexpected if fn is NULL)
    }

    # etas is initialized with length K-1, replicate(K - 1,  do.call(am$finit, ...))
    # Compared to the case above where etas is a matrix, here
    # there is no issue about how many 0-columns must be considered to ensure that
    # softmax_rows() returns K columns.
    #
    # Assume proper input from em_init()
    # stopifnot(length(etas) == K - 1)

    R <- length(etas) + 1L
    Maux <- matrix(0.0, nrow(features), R)

    for (r in seq_len(R - 1))
      Maux[,r+1] <- do.call(fn,
        c(list(X = features, etas = etas[[r]]), args))$Y

    rprobs <- softmax_rows(Maux)

  } else
    stop(sprintf("Unexpected class of argument 'etas': %s.", class(etas)))

  if (debug)
  {
    if (!all(abs(rowSums(rprobs) - 1) < .Machine$double.eps^0.8))
      stop("Some of the responsibilities do not sum to one.")
    if (any(rprobs < -1e-12)) stop("Some responsibilities turned out negative.")
    #if (nrow(rprobs) != N) stop("Wrong number of rows in responsibilities.")
    #if (ncol(rprobs) != K) stop("Wrong number of columns in responsibilities.")
  }

  rprobs
}

pretrain_amortization_kmeans <- function(model,
  nstart = 20, prob = 0.95,
  by_column = 1L,
  weight = 1e-4, # penalty weight
  optim_args = list(method = "BFGS", control = list(maxit = 200)),
  full_output = TRUE,
  debug = FALSE)
{
  if (!model$is.amortized)
    stop("Amortization features are not defined in 'model'.")

  if (!is.list(model$etas))
    stop("This function currently assumes neural-network parameters ",
         "stored as a list.")

  X <- model$X
  N <- model$N
  K <- model$K

  #if (K <= 1L)
  #  stop("'K' must be greater than one.")
  if (prob <= 1 / K || prob >= 1)
    stop("'prob' must lie between 1/K and 1.")

  # Flatten initial network parameters;
  # optim() takes as input a vector of parameters.
  # Get the structure of the list of parameters to
  # build model$etas within objective() defined below
  # and to build the output.

  etas_struct <- get_pars_structure_v2(model$etas)
  vec_etas_init <- unlist(model$etas, use.names = FALSE)

  # Initialize AR coefficients.
  #
  # Even if the model is not dynamic,
  # 'B_init' and 'E' need to be initialized.

  tmp <- em_init_arcoefs(model)
  B_init <- tmp$B
  E <- tmp$E

  # K-means

  km <- try(kmeans(E, centers = K, nstart = nstart), silent = TRUE)

  if (inherits(km, "try-error")) {
    #warning("K-means failed. Random initialization of responsibilities was applied.")
    warning("K-means failed. Try other options passed to pretrain_amortization_kmeans().")
  }

  # Give the components a reproducible order according to the
  # selected coordinate of the k-means centers.

  component_order <- order(km$centers[,by_column], decreasing = FALSE)
  cluster <- match(km$cluster, component_order)

  # Construct smoothed target responsibilities.
  # These are the probabilities to be matched with
  # the output of the network (after training it via optim()).

  target_rprobs <- matrix((1 - prob) / (K - 1), N, K)
  target_rprobs[cbind(seq_len(N), cluster)] <- prob

  # Cross-entropy objective.

  objective <- function(vec_etas)
  {
    eps <- 1e-12
    etas <- vec_to_list_v2(vec_etas, etas_struct)

    rprobs <- rprobs_from_amortization(etas = etas, K = K,
      amortization = model$amortization, debug = FALSE)

    cross_entropy <- -sum(target_rprobs * log(pmax(rprobs, eps))) / N
    penalty <- 0.5 * weight * mean(vec_etas^2)

    cross_entropy + penalty
  }

  # Fit the network to the k-means responsibilities.

  opt <- do.call(stats::optim,
    c(list(par = vec_etas_init, fn = objective), optim_args))

  etas <- vec_to_list_v2(opt$par, etas_struct)

  if (!full_output)
  {
    # when called from em_init_resp(),
    # rprobs_from_amortization() is done there
    return(etas)
  }

  # Compute responsibilities inferred by the amortization
  # function with the estimated parameters 'etas'.

  rprobs <- rprobs_from_amortization(etas = etas, K = K,
    amortization = model$amortization,
    debug = debug)

  # Assess how close the assignments learned by the
  # network match those obtained by k-means.

  agreement <- mean(max.col(rprobs) == cluster)

  list(etas = etas, rprobs = rprobs,
    target_rprobs = target_rprobs,
    cluster = cluster, kmeans = km,
    optim = opt, agreement = agreement)
}
