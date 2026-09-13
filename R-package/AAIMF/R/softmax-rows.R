
softmax_rows <- function(A)
{
  # Row-wise softmax: Amax is sweeped from A for numerical stability
  # See 'safe softmax method' in
  # https://en.wikipedia.org/wiki/Softmax_function#Numerical_algorithms
  Amax <- apply(A, 1L, max)
  exp_Acentered <- exp(sweep(A, 1, Amax, "-"))
  sweep(exp_Acentered, 1L, rowSums(exp_Acentered), "/")
}
