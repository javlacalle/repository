#
# Auxiliary functions to summarize clustering results.
#

summarize_component_sizes <- function(fit)
{
  rprobs <- fit$model$rprobs
  z_hat <- max.col(rprobs)

  data.frame(
    cluster = seq_len(ncol(rprobs)),
    posterior_pi = round(fit$posteriors$pis, 3),
    count = as.integer(tabulate(z_hat, nbins = ncol(rprobs))),
    effective_count = round(colSums(rprobs), 1))
}

cluster_mean_histograms <- function(X_hist, rprobs)
{
  K <- ncol(rprobs)
  out <- matrix(NA, nrow = K, ncol = ncol(X_hist))

  for (k in seq_len(K)) {
    w <- rprobs[,k]
    out[k, ] <- colSums(X_hist * w) / sum(w)
  }

  rownames(out) <- paste0("cluster_", seq_len(K))
  out
}

plot_cluster_rgb_hist <- function(h, bins_per_channel = 32L,
  txt = "Average RGB histogram") #main = "Average RGB histogram"
{
  idx_R <- seq_len(bins_per_channel)
  idx_G <- bins_per_channel + seq_len(bins_per_channel)
  idx_B <- 2L * bins_per_channel + seq_len(bins_per_channel)

  ylim <- range(h)

  plot(idx_R, h[idx_R], type = "n", ylim = ylim, col = "red",
       xlab = "Histogram bin", ylab = "Average proportion", main = "")
  mtext(txt, side=3, adj=0, cex=2)
  grid()

  lines(idx_R, h[idx_R], col = "red") #lty = 1
  lines(idx_R, h[idx_G], col = "green") #lty = 2
  lines(idx_R, h[idx_B], col = "blue") #lty = 3

  #legend("topright", legend = c("R", "G", "B"), lty = c(1, 2, 3), bty = "n")
}
