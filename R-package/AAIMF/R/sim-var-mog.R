var_mog_sim <- function(n,
  model = list(
    A = list(
      matrix(c(0.55, 0.10,
               0.05, 0.35), nrow = 2, byrow = TRUE),
      matrix(c(-0.25, 0.00,
                0.00, -0.10), nrow = 2, byrow = TRUE)),
    means = rbind(c(0, 0),c(3, 3)),
    Sigmas = list(
      diag(c(1/4, 1/4)),
      diag(c(1, 1))),
    probs = c(0.7, 0.3)),
  n.start = NA,
  seed = NULL)
{
  # NOTE documentation:
  # The first simulated values are initialized at zero and a burn-in period is discarded.
  # For highly persistent VAR models, it may be useful to increase n.start.

  # ------------------------------------------------------------
  # Model checks.
  # ------------------------------------------------------------

  if (is.null(model$A) || !is.list(model$A) || length(model$A) == 0L)
    stop("'model$A' must be a non-empty list of VAR coefficient matrices.")

  p <- length(model$A)
  D <- nrow(model$A[[1L]])

  if (D != ncol(model$A[[1L]]))
    stop("Each matrix in 'model$A' must be square.")

  for (i in seq_len(p))
  {
    if (!is.matrix(model$A[[i]]))
      stop("Each element of 'model$A' must be a matrix.")

    if (!all(dim(model$A[[i]]) == c(D, D)))
      stop("All matrices in 'model$A' must have dimension D x D.")
  }

  if (is.null(model$means))
    stop("'model$means' is missing or NULL.")

  if (is.vector(model$means))
    model$means <- matrix(model$means, nrow = 1L)

  if (!is.matrix(model$means))
    stop("'model$means' must be a matrix with one row per mixture component.")

  K <- nrow(model$means)

  if (ncol(model$means) != D)
    stop("The number of columns of 'model$means' must match the VAR dimension D.")

  if (is.null(model$probs))
    stop("'model$probs' is missing or NULL.")

  if (length(model$probs) != K)
    stop("The length of 'model$probs' must equal the number of mixture components.")

  if (any(model$probs < 0))
    stop("'model$probs' cannot contain negative values.")

  if (abs(sum(model$probs) - 1) > 1e-8)
    stop("The sum of 'model$probs' must be 1.")

  if (n <= 0L)
    stop("'n' must be strictly positive.")

  if (is.null(seed))
    stop("A non-NULL value for 'seed' was not set.")

  if (!is.null(seed))
    set.seed(seed)

  # ------------------------------------------------------------
  # Covariance / precision handling.
  # ------------------------------------------------------------

  if (!is.null(model$Sigmas))
  {
    if (!is.list(model$Sigmas) || length(model$Sigmas) != K)
      stop("'model$Sigmas' must be a list with one covariance matrix per component.")

    Sigmas <- model$Sigmas

    for (k in seq_len(K))
    {
      if (!is.matrix(Sigmas[[k]]) || !all(dim(Sigmas[[k]]) == c(D, D)))
        stop("Each element of 'model$Sigmas' must be a D x D matrix.")

      if (!isSymmetric(Sigmas[[k]]))
        stop("Each covariance matrix in 'model$Sigmas' must be symmetric.")

      eigvals <- eigen(Sigmas[[k]], symmetric = TRUE, only.values = TRUE)$values

      if (min(eigvals) <= 0)
        stop("Each covariance matrix in 'model$Sigmas' must be positive definite.")
    }
  } else if (!is.null(model$Lambdas)) {
    if (!is.list(model$Lambdas) || length(model$Lambdas) != K)
      stop("'model$Lambdas' must be a list with one precision matrix per component.")

    Sigmas <- vector("list", K)

    for (k in seq_len(K))
    {
      if (!is.matrix(model$Lambdas[[k]]) || !all(dim(model$Lambdas[[k]]) == c(D, D)))
        stop("Each element of 'model$Lambdas' must be a D x D matrix.")

      if (!isSymmetric(model$Lambdas[[k]]))
        stop("Each precision matrix in 'model$Lambdas' must be symmetric.")

      eigvals <- eigen(model$Lambdas[[k]], symmetric = TRUE, only.values = TRUE)$values

      if (min(eigvals) <= 0)
        stop("Each precision matrix in 'model$Lambdas' must be positive definite.")

      Sigmas[[k]] <- solve(model$Lambdas[[k]])
    }
  } else {
    stop("Either 'model$Sigmas' or 'model$Lambdas' must be provided.")
  }

  # ------------------------------------------------------------
  # Stationarity check through the VAR companion matrix.
  # The VAR is stable if all eigenvalues of the companion matrix
  # have modulus strictly less than one.
  # ------------------------------------------------------------

  companion <- matrix(0, D * p, D * p)

  companion[seq_len(D), ] <- do.call(cbind, model$A)

  if (p > 1L) {
    companion[(D + 1L):(D * p), seq_len(D * (p - 1L))] <- diag(D * (p - 1L))
  }

  max_root <- max(Mod(eigen(companion, only.values = TRUE)$values))

  if (max_root >= 1)
    stop("VAR model is not stationary: the companion matrix has eigenvalue modulus >= 1.")

  # ------------------------------------------------------------
  # Burn-in length.
  # ------------------------------------------------------------

  if (is.na(n.start))
    n.start <- max(100L, 10L * p)

  if (n.start < p)
    stop("Burn-in 'n.start' must be at least p.")

  n.total <- n.start + n

  # ------------------------------------------------------------
  # Simulate mixture states and innovations.
  # ------------------------------------------------------------

  states <- sample(K, size = n.total, replace = TRUE, prob = model$probs)

  eps <- matrix(NA_real_, n.total, D)

  for (t in seq_len(n.total))
  {
    k <- states[t]
    #eps[t,] <- as.numeric(MASS::mvrnorm(n = 1L, mu = model$means[k, ], Sigma = Sigmas[[k]]))
    eps[t,] <- as.numeric(mvrnorm(n = 1L, mu = model$means[k, ], Sigma = Sigmas[[k]]))
  }

  # ------------------------------------------------------------
  # Recursive VAR simulation.
  # Initial values are zero. They are discarded by burn-in.
  # ------------------------------------------------------------

  Xfull <- matrix(0, nrow = n.total, ncol = D)

  for (t in seq_len(n.total))
  {
    xhat <- rep(0, D)

    for (i in seq_len(p))
    {
      if (t - i > 0L)
        xhat <- xhat + model$A[[i]] %*% Xfull[t - i, ]
    }

    Xfull[t, ] <- as.numeric(xhat + eps[t, ])
  }

  # ------------------------------------------------------------
  # Remove burn-in.
  # ------------------------------------------------------------

  X <- Xfull[(n.start + 1L):n.total, , drop = FALSE]
  states <- states[(n.start + 1L):n.total]

  # ------------------------------------------------------------
  # Build lag matrix for later EM implementation.
  # This assumes that 'build_var_lags_matrix()' returns:
  #   Y:     n - p by D response matrix
  #   H:     n - p by Dp lag matrix
  #   Xlags: same as H, if that name is preferred elsewhere
  # ------------------------------------------------------------

  tmp <- build_var_lags_matrix(X, p)

  c(mget(names(formals(var_mog_sim)), envir=environment()),
  list(
    X = tmp$Y,
    Xlags = tmp$Xlags,
    H = tmp$Xlags,
    raw_X = X,
    states = states[(p + 1L):n],
    raw_states = as.ts(states),
    #model = model,
    #Sigmas = Sigmas,
    K = K,
    D = D,
    p = p,
    companion_roots = eigen(companion, only.values = TRUE)$values
  ))
}

if (interactive())
{
dgp <- var_mog_sim(
  n = 500,
  model = list(
    A = list(
      matrix(c(0.55, 0.10,
               0.05, 0.35), nrow = 2, byrow = TRUE),
      matrix(c(-0.25, 0.00,
                0.00, -0.10), nrow = 2, byrow = TRUE)
    ),
    means = rbind(
      c(0, 0),
      c(3, 3)
    ),
    Sigmas = list(
      diag(c(1 / 4, 1 / 4)),
      diag(c(1, 1))
    ),
    probs = c(0.7, 0.3)
  ),
  seed = 123
)

dim(dgp$X)
dim(dgp$Xlags)
table(dgp$states)
dgp$D
dgp$p
dgp$K
}
