
make_model <- function(
  family = c('mog_normal_wishart', 'mog_normal_gamma', 'var_mog_normal_wishart', 'var_mog_normal_gamma',
             'mog_full_cov', 'mog_diagonal_cov', 'var_mog_full_cov', 'var_mog_diagonal_cov'),
  X, K,
  p = 0L, Xlags = NULL, # var_mog_sim() returns Xlags; in simulations use it avoid repeat; in practice use p
  priors = list(),
  amortization = list(features_fn = NULL, fn = NULL, args_fn = list()), nn_width=4)
  #init_amortization = list(finit = NULL, init_args = list())
{
  # TODO see define args_fn within args_fn()

  # Creates an object of the class 'var_mog' containing the
  # information required by the package's functions.
  #
  # Using this class standardizes argument passing and avoids
  # redundant validation checks inside numerically intensive functions.

  # NOTE Documentation:
  # make_model() sets the priors parameters;
  # mog_init() initializes the responsibilities and the model parameters.

  # NOTE Documentation:
  # The prototype of 'features_fn' must be function(X), where X will be taken
  # to be the argument 'X'. Example: features_fn = function(X) cbind(1, X, X^2).
  # In the future further arguments to this function may be allowed.
  #
  # It was convenient to pass the function 'features_fn' instead of
  # passing the matrix of features themselves. In this way, the model object
  # returned by make_model() keeps the function that generates the features
  # and can be used for example by fitted.avem_fit().
  #
  # The model object includes in the slot 'amortization' an element named 'features'
  # containing features_fn(X); so the features themselves are also stored in the
  # object and need not be recomputed each time a function needs them.
  #
  # WARNING: When calling these elements, be aware that:
  # model$amortization$features contains the matrix of features, not features_fn;
  # model$amortization$feat returns NULL, as the name of 2 elements start by 'feat'.

  family <- match.arg(family)

  if (!is.numeric(K) || K %% 1 != 0 || K <= 1)
    stop("'K' must be an integer greater than 1.")

  # X

  if (!is.matrix(X))
    stop("'X' must be a matrix.")

  N <- nrow(X)
  D <- ncol(X)

  if (N < K)
    stop(sprintf("The number of observations is lower than the number of components, N = %d < K = %d.", N, K))

  # Amortization.

  feats_fn <- amortization$features_fn
  is.amortized <- !is.null(feats_fn)

  if (is.amortized)
  {
    feats <- feats_fn(X)
    amortization$features <- feats

    if (!is.matrix(feats) || nrow(feats) != N)
      stop("'amortization$features_fn' must be a function returning a matrix with 'nrow(X)' rows.")

    if (amortization$fn == "forward_network")
    {
      if (is.null(amortization$nn_width))
      {
        # NOTE documentation:
        # Do not set 'nn_with' within 'amortization$args_fn';
        # those are the arguments passed to the inference mapping,
        # at that point the design of the neural-network is already
        # set in the parameters of the network 'model$etas'.

        # The argument 'nn_width' is a bit out of place,
        # it is a parameter required to initialize the
        # neural-network, not to define the amortization
        # features or the inference function.
        # By default init_network_random() uses width = 4,
        # so this could be skipped.
        #
        # Eventually it has been decided to keep 'nn_width' there
        # here because it is a relevant parameter. It helps to
        # make sure that we are aware of it when defining the model
        # (either because it is explicitly passed to make_model()
        # or because it is reminded by the warning below).
        #
        # In the future, an additional argument 'init_amortization'
        # as used in a previous version could be added. For now,
        # this argument can be avoided, simplifying the usage
        # of this function. Just keep 'nn_width' for the reason above.

        warning("The shape of the neural-network was not specified. ",
                "amortization$nn_width = 4 was set.")
        amortization$nn_width <- 4
      }
    }
  }

  # Xlags

  if (!is.numeric(p) || p %% 1 != 0 || p < 0) # p != as.integer(p)
    stop("'p' must be a positive integer or zero.")

  if (p > 0)
  {
    #if (family %in% c('mog_full_cov', 'mog_diagonal_cov'))
    #  stop(sprintf("'p' > 0 with family '%s% is not implemented in this interface. Use 'var_mog_em' instead.", family))

    if (is.null(Xlags))
    {
      tmp <- if (ncol(X) == 1) build_lags_matrix(X, p) else build_var_lags_matrix(X, p)
      Xlags <- tmp$Xlags
      # enforce nrow(X) == nrow(Xlags)
      X <- tmp$X
      N <- nrow(X)
      if (is.amortized)
        amortization$features <- feats[-seq_len(p),,drop=FALSE]

    } else {
     if (nrow(X) != nrow(Xlags))
       stop("'X' and 'Xlags' must have the same number of rows. ",
           "Missing observations due to lags should be excluded.")

     if (ncol(Xlags) != D * p)
       stop("'Xlags' must contain p x ncol(X) columns ('p' lags for each variable in 'X').")
    }
  }

  # Priors: Complete (if required) with default prior parameters.

  priors_full <- switch(family,
    "var_mog_normal_wishart" = list(lambda0 = 1, alpha_a0 = 0.01, alpha_r0 = 0.01,
          m0 = 0, kappa0 = 0.01, W0 = diag(D), nu0 = D + 2),

    "var_mog_normal_gamma" = list(lambda0  = 1, alpha_a0 = 0.01, alpha_r0 = 0.01,
          m0 = 0, kappa0 = 0.01, a0 = 1, b0 = 1),

    "mog_normal_wishart" = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL),

    "mog_normal_gamma" = list(alpha0 = 1, beta0 = 1, m0 = NULL, a0 = 1, b0 = 1),

    list() # default, family = 'mog_full_cov', 'mog_diagonal_cov', 'var_mog_full_cov', 'var_mog_diagonal_cov'
  )

  if (length(priors))
  {
    if (family %in% c('mog_full_cov', 'mog_diagonal_cov', 'var_mog_full_cov', 'var_mog_diagonal_cov'))
      stop(sprintf("Argument 'priors' should be left empty. It does not apply to family '%s'.", family,))

    if (!all(names(priors) %in% names(priors_full)))
        stop(sprintf("The following prior parameters do not belong to '%s':\n    %s.\nUse any of\n    %s.",
            family, paste(setdiff(names(priors), names(priors_full)), collapse=", "),
            paste(names(priors_full), collapse=", ")))

    priors_full[names(priors)] <- priors
  }

  # Output.

  is.bayesian <- family %in% c('mog_normal_wishart', 'mog_normal_gamma', 'var_mog_normal_wishart', 'var_mog_normal_gamma')
  is.fullcov <- grepl("full_cov", family) | grepl("normal_wishart", family)

  model <- list(family = family,
    X = X, Xlags = Xlags,
    #priors = priors_full,
    etas = NULL, rprobs = NULL, pars = NULL,
    amortization = amortization, #init_amortization = init_amortization,
    K = K, N = N, D = D, p = p,
    is.classical = !is.bayesian,
    is.bayesian = is.bayesian,
    is.dynamic = p > 0,
    is.diagonal = !is.fullcov,
    is.fullcov = is.fullcov,
    is.amortized = is.amortized)

  if (length(priors_full) > 0)
    model <- append(model, list(priors = priors_full), after = match("Xlags", names(model)))

  # TODO print() method for this class
  class(model) <- "var_mog"
  model
}
