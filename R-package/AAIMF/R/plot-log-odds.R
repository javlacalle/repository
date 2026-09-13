##
## Compute fitted values for the log-odds.
##
## NOTE
## The fitted values and the plots
## are currently implemented only for the 'mog_diagonal_cov' model
## with K = 2, D = 2.
##


# FIXME add label for type!="difference" for 0-level contour: estimated classification boundary
#
# Notebook, interpretation
#The exact and amortized boundaries both lie in the gap between the two observed groups. The amortized boundary is more nearly #vertical, whereas the exact boundary is oblique. The discrepancy is more pronounced away from the observations.
#
#This suggests that the network learned a useful classifier but did not reproduce the exact log-odds surface globally. It does #not look anomalous for a numerically optimized amortization map.
#
#The exact and amortized boundaries intersect near the sparsely populated region between the groups but have different curvature. Both still separate most observations correctly.
#This is consistent with a good classification mapping but a non-negligible inference gap. Since the ELBO depends on the responsibility probabilities rather than only on hard classifications, this difference can matter even when accuracy remains high.
#
#The red and blue boundaries are relatively close around the region separating the observations. They diverge more strongly in regions with little or no data. This is a fairly natural outcome: the network is fitted through the observed feature values and is not constrained to reproduce the analytical posterior surface far outside the data cloud.
#
#The white equality contour crossing other parts of the plane does not indicate an error. It merely locates points where the two log-odds happen to agree.
#
#Iteration 33 This plot is qualitatively different. Neither the exact nor the amortized zero boundary appears inside the displayed region, and the difference has essentially one sign throughout the grid.

#object = fit2_elbo_ir
fitted.avem_fit <- function(object, ...,
  ngrid = 200L, fgrid = 0.5)
  #features_fn = function(X) cbind(1, X, X^2))
{
  #fit <- object
  if (object$model$family != "mog_diagonal_cov")
    stop("The method 'fitted' is currently implemented only for family 'mog_diagonal_cov'.")

  if (object$model$K != 2)
    stop("The method 'fitted' is currently implemented only for K=2 componentes.")

  if (object$model$D != 2)
    stop("The method 'fitted' is currently implemented only for D=2 variables.")

  model <- object$model
  X <- model$X
  pars <- model$pars
  features_fn <- object$model$amortization$features_fn
  is_amortized <- object$model$is.amortized

  # Compute the coefficients of the analytical expression of the log-odds
  # for the 'mog_diagonal_cov' with K = 2, D = 2.

  lambda1 <- as.numeric(pars$Lambdas[[1]])
  lambda2 <- as.numeric(pars$Lambdas[[2]])

  mu1 <- pars$mus[1,]
  mu2 <- pars$mus[2,]

  # coefficients of the theoretical log-odds ('mog_diagonal_cov' model)
  c0 <- log(pars$pis[2] / pars$pis[1]) +
    0.5 * sum(log(lambda2 / lambda1)) -
    0.5 * sum(lambda2 * mu2^2 - lambda1 * mu1^2)
  a <- lambda2 * mu2 - lambda1 * mu1
  b <- 0.5 * (lambda1 - lambda2)

  coefs <- list(c = c0, a = a, b = b)

  # Compute the fitted log-odds for the observed data
  # according to the exact analytical expression and
  # the estimated parameter values (observed_exact).

  observed_exact <- as.numeric(c0 + X %*% a + X^2 %*% b)
  obs_prob_2 <- plogis(observed_exact) # second component

  if (ngrid == 0)
    return(list(coefs = coefs, log_odds = observed_exact,
      obs_prob_2 = obs_prob_2))

  # from this point ngrid is > 0 (no need for an if-block)

  #x1 <- seq(min(X[,1]) - 0.5, max(X[,1]) + 0.5, length.out = ngrid)
  #x2 <- seq(min(X[,2]) - 0.5, max(X[,2]) + 0.5, length.out = ngrid)
  Xmin <- apply(X, 2L, min)
  Xmax <- apply(X, 2L, max)
  x1 <- seq(Xmin[1] - fgrid * abs(Xmin[1]), Xmax[1] + fgrid * abs(Xmax[1]), length.out = ngrid)
  x2 <- seq(Xmin[2] - fgrid * abs(Xmin[2]), Xmax[2] + fgrid * abs(Xmax[2]), length.out = ngrid)
  Grid <- as.matrix(expand.grid(X1 = x1, X2 = x2))

  # Compute the exact fitted log-odds for a grid of values

  grid_exact <- as.numeric(coefs$c + Grid %*% coefs$a + Grid^2 %*% coefs$b)
  Z_exact <- matrix(grid_exact, length(x1), length(x2))

  if (is_amortized)
  {
    etas <- object$model$etas[[1]]

    features_observed <- object$model$amortization$features #features_fn(X)
    observed_amortized <- as.numeric(forward_network(features_observed, etas)$Y)

    features_grid <- features_fn(Grid)
    grid_amortized <- as.numeric(forward_network(features_grid, etas)$Y)
    Z_amortized <- matrix(grid_amortized, length(x1), length(x2))

    observed_difference <- observed_amortized - observed_exact
    grid_difference <- grid_amortized - grid_exact
    Z_difference <- matrix(grid_difference, length(x1), length(x2))

  } else {
    observed_amortized <- NULL
    Z_amortized <- NULL
    observed_difference <- NULL
    Z_difference <- NULL
  }

  res <- list(x1 = x1, x2 = x2, X = X,
    coefs = coefs,
    observed_exact = observed_exact, observed_amortized = observed_amortized,
    observed_difference = observed_difference,
    grid_exact = Z_exact, grid_amortized = Z_amortized,
    grid_difference = Z_difference)

  class(res) <- "fitted_log_odds"
  res
}

plot_fitted_log_odds <- function(x,
  type = c("exact", "amortized", "difference"),
  main = NULL, pcols = NULL)
{
  # NOTE for future
  # see color data points depending on rprobs assigned to each of them

  stopifnot(class(x) == "fitted_log_odds")

  type <- match.arg(type)
  is_amortized <- !is.null(x$observed_amortized)

  if (type == "exact")
  {
    Z <- x$grid_exact
    observed_values <- x$observed_exact

    if (is_amortized) {
      zlim <- range(c(x$grid_exact, x$grid_amortized), finite = TRUE)
    } else
      zlim <- range(x$grid_exact, finite = TRUE)

  } else
  if (type == "amortized")
  {
    Z <- x$grid_amortized
    observed_values <- x$observed_amortized

    zlim <- range(c(x$grid_exact, x$grid_amortized), finite = TRUE)

  } else
  if (type == "difference")
  {
    Z <- x$grid_difference
    observed_values <- x$observed_difference

    diff_lim <- max(abs(Z), na.rm = TRUE)
    #if (diff_lim == 0)
    #  diff_lim <- .Machine$double.eps
    zlim <- c(-diff_lim, diff_lim)
  }

  if (is.null(main)) {
    main = switch(type,
      "exact" = "Exact log-odds based on parameter estimates",
      "amortized" = "Amortized log-odds",
      "difference" = "Amortized minus exact log-odds")
  }

  palette <- hcl.colors(101, palette = "Cyan-Magenta") #"Blue-Red 3"

  # If type != "difference", the following displays the exact/amortize log-odds.
  # If type = "difference", the following displays the signed difference between the log-odds.

  image(x$x1, x$x2, Z, col = palette, zlim = zlim,
    xlab = expression(x[1]), ylab = expression(x[2]),
    #main = sprintf("Exact posterior log-odds based on parameter estimates. Iter =  %s", iter))
    main = main)

  # Data points (the sample).

  if (is.null(pcols)) {
    points(x$X[,1], x$X[,2], pch = 16, col = "black", cex = 1.2)
  } else {
    stopifnot(length(pcols) == nrow(x$X))
    points(x$X[,1], x$X[,2], pch = 16, col = pcols, cex = 1.2)
  }

  # Contours.

  if (type != "difference")
  {
    # 0-level contour: stands for a responsibility probability for component 2 equal to 1/2.
    contour(x$x1, x$x2, Z, levels = 0, add = TRUE, drawlabels = FALSE, lwd = 3, col = "black")

    legend("topleft",
      legend = c("Odds favour component 1", "Odds favour component 2", "Classification boundary p = 1/2"),
      fill = c(palette[1], tail(palette, 1), NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
      border = c(NA, NA, NA), bty = "n")

  } else { # type == "difference"
    # 0-level contour: the exact and amortized log-odds agree.
    contour(x$x1, x$x2, Z, levels = 0, add = TRUE, drawlabels = FALSE, lwd = 3, col = "white")

    # Classification boundaries.
    contour(x$x1, x$x2, x$grid_exact, levels = 0, add = TRUE, drawlabels = FALSE, lwd = 3, col = "red")
    contour(x$x1, x$x2, x$grid_amortized, levels = 0, add = TRUE, drawlabels = FALSE, lwd = 3, col = "blue")

    legend("topleft",
      legend = c("Exact boundary p = 1/2", "Amortized boundary p = 1/2", "Equal log-odds"),
      lty = c(1, 1, 1), lwd = c(3, 3, 2),
      col = c("red", "blue", "white"), bty = "n")
  }

  #invisible(res)
}
