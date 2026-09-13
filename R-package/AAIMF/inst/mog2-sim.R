library(AAIMF)

set.seed(123)

niter <- 20
N <- 500
K <- 2

res <- vector("list", niter)

for (iter in seq_len(niter))
{
  dgp <- mog2_sim(N = N,
    pis = c(0.35, 0.65),
    mus = rbind(-2, 2),
    Sigmas = list(matrix(0.4, 1, 1), matrix(0.8, 1, 1)),
    seed = 1000 + iter)

  ## Pure CAVI.
  t_cavi <- system.time({
    fit_cavi <- mog2_cavi(X = dgp$X, K = dgp$K,
      init_args = list(seed = 2000 + iter),
      maxiter = 500, tol = 1e-8,
      verbose = FALSE, debug = FALSE)
  })

  ## Linear amortization.
  am_feats <- cbind(1, dgp$X)
  #am_feats <- cbind(1, dgp$X, dgp$X^2)

  init_lin <- mog2_init(X = dgp$X, K = dgp$K,
    fixed_pars = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL),
    amortization = list(features = am_feats, init_args = list(sd = 0.01)),
    seed = 3000 + iter)

  model_struct_lin <- get_pars_structure(init_lin$pars)

  vec_etas_lin <- if (is.list(init_lin$pars$etas)) {
    unlist(init_lin$pars$etas, use.names = FALSE)
  } else {
    as.vector(init_lin$pars$etas)
  }

  t_lin <- system.time({
    fit_lin <- mog2_cavi_amortized(vec_etas = vec_etas_lin,
      model_struct = model_struct_lin, X = dgp$X,
      fixed_pars = init_lin$fixed_pars,
      amortization = list(features = am_feats, fn = function(X, etas) X %*% etas),
      optim_args = list(method = "BFGS", hessian = FALSE, control = list(maxit = 1000)))
  })

  ## Neural-network amortization.
  am_feats_nn <- cbind(dgp$X)

  init_nn <- mog2_init(X = dgp$X, K = dgp$K,
    fixed_pars = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL),
    amortization = list(features = am_feats_nn,
      fn = 'forward_network', args_fn = list(return_cache = FALSE),
      finit = 'init_network', init_args = list(input_dim = ncol(am_feats_nn),
        width = 10, output_dim = dgp$D, seed = 4000 + iter)))

  model_struct_nn <- get_pars_structure(init_nn$pars)

  vec_etas_nn <- unlist(init_nn$pars$etas, use.names = FALSE)

  t_nn <- system.time({
    fit_nn <- mog2_cavi_amortized(vec_etas = vec_etas_nn, model_struct = model_struct_nn,
      X = dgp$X, fixed_pars = init_nn$fixed_pars,
      amortization = list(features = am_feats_nn,
        fn = forward_network, args = list(return_cache = FALSE)),
      optim_args = list(method = "BFGS", hessian = FALSE, control = list(maxit = 1000)))
  })

  cavi_elbo <- tail(fit_cavi$elbo_path, 1)

  res[[iter]] <- data.frame(iter = iter,
    method = c("CAVI", "Linear AVI", "NN AVI"),
    elbo = c(cavi_elbo, fit_lin$elbo, fit_nn$elbo),
    amortization_gap = c(0, cavi_elbo - fit_lin$elbo, cavi_elbo - fit_nn$elbo),
    time_sec = c(unname(t_cavi["elapsed"]), unname(t_lin["elapsed"]), unname(t_nn["elapsed"])),
    mu_1 = c(sort(fit_cavi$pars$mu_ms[,1])[1], sort(fit_lin$pars$mu_ms[,1])[1], sort(fit_nn$pars$mu_ms[,1])[1]),
    mu_2 = c(sort(fit_cavi$pars$mu_ms[,1])[2], sort(fit_lin$pars$mu_ms[,1])[2], sort(fit_nn$pars$mu_ms[,1])[2]))
}

res_df <- do.call(rbind, res)

aggregate(
  cbind(elbo, amortization_gap, time_sec, mu_1, mu_2) ~ method,
  data = res_df,
  FUN = function(x) c(mean = mean(x), sd = sd(x)))
