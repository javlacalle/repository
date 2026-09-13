
# ----------------------------------------------------------------
# Neural Network.
# Amortization mapping.
# ----------------------------------------------------------------

init_network_random <- function(input_dim = 1, width = 4, output_dim = 1, seed = NULL)
{
  if (!is.null(seed))
    set.seed(seed)

  list(
    W1 = matrix(rnorm(input_dim * width, sd = sqrt(2 / input_dim)), input_dim, width),
    b1 = matrix(0, 1, width),

    W2 = matrix(rnorm(width * width, sd = sqrt(2 / width)), width, width),
    b2 = matrix(0, 1, width),

    W3 = matrix(rnorm(width * output_dim, sd = sqrt(2 / width)), width, output_dim),
    b3 = matrix(0, 1, output_dim)
  )
}

forward_network <- function(X, etas, return_cache = FALSE)
{
  # X: n x 1
  Z1 <- X %*% etas$W1 + matrix(etas$b1, nrow = NROW(X), ncol = NCOL(etas$b1), byrow = TRUE)
  # relu
  #FIXME see safer version for matrices
  A1 <- pmax(Z1, 0)

  Z2 <- A1 %*% etas$W2 + matrix(etas$b2, nrow = NROW(X), ncol = NCOL(etas$b2), byrow = TRUE)
  # relu
  A2 <- pmax(Z2, 0)

  Y <- A2 %*% etas$W3 + matrix(etas$b3, nrow = NROW(X), ncol = NCOL(etas$b3), byrow = TRUE)

  cache <- if (return_cache) list(X = X, Z1 = Z1, A1 = A1, Z2 = Z2, A2 = A2) else NULL

  list(Y = Y, cache = cache)
}

# for vectors
#relu_grad <- function(z) {
#  as.numeric(z > 0)
#}

relu_grad <- function(z) {
  # this version is safer for matrices
  res <- z
  res[,] <- as.numeric(z > 0)
  res
}

backward_network <- function(net, cache, dY)
{
  # dY is the derivative of ELBO with respect to network_output (N x 1).
  # Returns gradients with the same structure as net.

  X <- cache$X
  Z1 <- cache$Z1
  A1 <- cache$A1
  Z2 <- cache$Z2
  A2 <- cache$A2

  dW3 <- t(A2) %*% dY
  db3 <- matrix(colSums(dY), 1, ncol(dY))

  dA2 <- dY %*% t(etas$W3)
  dZ2 <- dA2 * relu_grad(Z2)

  dW2 <- t(A1) %*% dZ2
  db2 <- matrix(colSums(dZ2), 1, ncol(dZ2))

  dA1 <- dZ2 %*% t(etas$W2)
  dZ1 <- dA1 * relu_grad(Z1)

  dW1 <- t(X) %*% dZ1
  db1 <- matrix(colSums(dZ1), 1, ncol(dZ1))

  list(W1 = dW1, b1 = db1, W2 = dW2, b2 = db2, W3 = dW3, b3 = db3)
}
