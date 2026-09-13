
# TODO merge into a single function and homogenize output format.

# Generate a matrix of lagged versions of the input time series (vector) 'x'.

build_lags_matrix <- function(x, p)
{
  #stopifnot(length(x) > p, p >= 0)
  #if (p == 0)
  #  return(list(X = matrix(x), Xlags = NULL))
  #stopifnot(length(x) > p, p >= 1)

  n <- length(x)
  Xlags <- matrix(NA, nrow = n - p, ncol = p)

  for (i in seq_len(p)) {
    Xlags[,i] <- x[(p + 1 - i):(n - i)]
  }

  # It seems stats::embed does not speed up the loop above;
  #stats::embed(x, dimension = p + 1)[,-1,drop=FALSE]

  colnames(Xlags) <- paste0("lag", seq_len(p))
  list(X = matrix(x[-seq_len(p)]), Xlags = Xlags)
}

build_var_lags_matrix <- function(X, p)
{
  if (!is.matrix(X))
    X <- as.matrix(X)

  n <- nrow(X)
  D <- ncol(X)

  if (p < 1L)
    stop("'p' must be at least 1.")

  if (n <= p)
    stop("The number of rows of 'X' must be larger than 'p'.")

  Z <- embed(X, dimension = p + 1L)

  Y <- Z[, seq_len(D), drop = FALSE]
  Xlags <- Z[, -(seq_len(D)), drop = FALSE]

  colnames(Y) <- paste0("x", seq_len(D))
  colnames(Xlags) <- unlist(
    lapply(seq_len(p), function(i) paste0("x", seq_len(D), "_lag", i)))

  list(Y = Y, X = Y, Xlags = Xlags,
    p = p, D = D, n_used = n - p)
}
