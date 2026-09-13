#
# Print methods for classes defined in the package.
#


#TODO print method for avem_fit (reuse print from above, just do header indicating target fn)


print.em_fit <- function(x, digits = 4, ...)
{
  if (!inherits(x, "em_fit"))
    stop("'x' must be an object of class 'em_fit'.")

  #family <- get(x$model$family)
  model <- x$model
  K <- model$K
  p <- model$p
  pars <- model$pars
  Lambdas <- pars$Lambdas
  colSums_rprobs <- x$colSums_rprobs

  cat("\n")
  if (p == 0) {
    cat("Expectation-Maximization fit for mixture-of-Gaussians model\n")
  } else # p > 0 (Assumed proper input)
    cat("Expectation-Maximization fit for VAR mixture-of-Gaussians-noise model\n")

  #isDiagonal <- grepl("diagonal", model$family)
  isFullcov <- model$is.fullcov
  if (isFullcov) {
    cat("Full covariance matrix\n")
  } else
    cat("Diagonal covariance matrix\n")
  cat("====================================================\n")

  cat("\nModel dimensions\n")
  cat("----------------------------------------------------------------\n")
  cat("Number of observations:", model$N, "\n")
  cat("Number of variables D:", model$D, "\n")
  cat("Number of mixture components K:", K, "\n")
  if (p > 0)
    cat("AR order p:", p, "\n")

  cat("\nConvergence\n")
  cat("----------------------------------------------------------------\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$niter, "\n")

  if (!is.null(x$loglik_path))
  {
    loglik_used <- x$loglik_path[is.finite(x$loglik_path)]

    if (length(loglik_used) > 0)
    {
      cat("Initial Log-Likelihood:", round(loglik_used[1], digits), "\n")
      cat("Final Log-Likelihood:", round(tail(loglik_used, 1), digits), "\n")
    }
  }

  cat("\nFitted parameters summary \n")
  cat("----------------------------------------------------------------\n")

  cat("\nMixture probabilities:\n")
  print(round(pars$pis, digits))

  cat("\nEffective component sizes:\n")
  print(round(colSums_rprobs, digits))

  cat("\nComponent/Innovation means, m_k:\n")
  print(round(pars$mus, digits))

  if (p > 0)
  {
    cat("\nVAR coefficient matrices A_1, ..., A_p:\n")
    A_hat <- B_to_A(pars$B, D = model$D, p = p)
    invisible(lapply(seq_along(A_hat), function(j) {
      cat("\nA_", j, ":\n", sep = "")
      print(round(A_hat[[j]], digits))
    }))
  }

  cat("\nComponent/Innovation precisions:\n")
  #invisible(lapply(pars$Lambdas, function(x) print(round(x, digits))))
  invisible(lapply(Lambdas, function(x) print(round(x, digits))))

  cat("\nComponent/Innovation variances:\n")
  # TODO see use ridge factor.
  if (inherits(Lambdas[[1]], "matrix")) {
    invisible(lapply(Lambdas, function(x) print(round(solve(x), digits))))
  } else
    invisible(lapply(Lambdas, function(x) print(round(1/x, digits))))

  cat("\n")
  invisible(x)
}


print.em_experimental_fit <- function(x, digits = 4, ...)
{
  if (!inherits(x, "em_experimental_fit"))
    stop("'x' must be an object of class 'em_experimental_fit'.")

  y <- x
  class(y) <- "em_fit"

  # FIXME make it clear in the output; documentation: optim code (not TRUE/FALSE).
  y$converged <- x$optim$convergence
  # FIXME (function and gradient counts, not explicit in the printed output)
  y$niter <- x$optim$counts
  y$loglik_path <- NULL
  y$colSums_rprobs <- colSums(y$model$rprobs)

  cat("\n")
  cat("Amortized (hybrid) Expectation-Maximization:\n")
  cat("====================================================\n")
  print(y)

  cat("Final logLik:", round(x$loglik, digits), "\n")

  cat("\nAmortization function parameters.\n")
  cat("----------------------------------------------------------------\n")

  if (is.matrix(x$model$etas)) {
    print(round(x$model$etas, digits))
  } else{
    n_am_pars <- length(unlist(x$model$etas, TRUE, FALSE))
    if (n_am_pars > 20) {
      cat("The number of amortization parameters is too large to be printed:", n_am_pars, ".\n")
      cat("Inspect 'x$model$etas'.\n")
    } else {
      #cat("Inspect 'x$model$etas'; contains", n_am_pars, "parameters. \n")
      invisible(lapply(x$model$etas, function(x) print(round(x, digits))))
    }
  }
}
