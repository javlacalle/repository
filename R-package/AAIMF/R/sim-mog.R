
mog_sim <- function(N = 500,
  pis = c(0.4, 0.6), mus = rbind(c(-2, 0), c(2, 0)), Sigmas = NULL,
  seed = NULL)
{
  if (!is.null(seed))
    set.seed(seed)

  #mus <- as.matrix(mus)
  D <- ncol(mus)
  K <- length(pis)

  if (is.null(Sigmas)) {
    Sigmas <- replicate(K, diag(D), simplify = FALSE)
  } else
  if (length(Sigmas) != K)
    stop("'Sigmas' must be a list of K covariance matrices.")

  zs <- sample.int(K, size = N, replace = TRUE, prob = pis)
  X <- matrix(NA, nrow = N, ncol = D)

  for (k in seq_len(K))
  {
    idx <- which(zs == k)
    nk <- length(idx)
    if (nk > 0L)
    {
      L <- chol(Sigmas[[k]])
      X[idx, ] <- matrix(rnorm(nk * D), nk, D) %*% L +
        matrix(mus[k, ], nk, D, byrow = TRUE)
    }
  }
  res <-  list(X = X, states = zs,
               pis = pis, mus = mus, Sigmas = Sigmas,
               N = N, D = D, K = K)
  class(res) <- "mog_sim"
  res
}

print.mog_sim <- function(x, ...)
{
  cat("Simulation from the mixture of Gaussians model\n")
  sprintf("N = %d; D = %d, K = %d", x$N, x$D, x$K)
  cat("\nMixture probabilities:\n")
  print(x$pis)
  cat("\nComponent means:\n")
  print(x$mus)
  cat("\nClass counts:\n")
  print(table(factor(x$zs, levels = seq_len(x$K))))
  invisible(x)
}
