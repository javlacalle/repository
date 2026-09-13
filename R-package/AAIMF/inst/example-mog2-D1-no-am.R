library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(-2, 2), Sigmas = list(matrix(0.4, 1, 1), matrix(0.8, 1, 1)),
  seed = 123)
#X <- dgp$X
#print(dgp)

#dput(dgp, file="test-dgp.txt")

# Initialize the parameters of the model.

# family = mog_normal_wishart
# K = dgp$K
# Xlags = NULL
# fixed_pars = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL)
# amortization = list(features = NULL, type = c("linear", "network"), args = list())
# kmeans_nstart = 10
# kmeans_prob = 0.95
# seed = 123321

init <- mog_init(mog_normal_wishart, dgp$X, dgp$K,
  fixed_pars = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL),
  amortization = list(features = NULL),
  kmeans_nstart = 10, kmeans_prob = 0.95,
  seed = 123321)

#init$pars

#dput(as_bishop_gmm_init_elbo(init), file="test-init.txt")

model_struct <- get_pars_structure(init$pars)

# Evaluate the ELBO at the initial parameter values.

# family = mog_normal_wishart
# X = dgp$X
# pars = init$pars
# fixed_pars = init$fixed_pars
# iNDK = init$iNDK
# amortization = list(features = NULL, fn = NULL, args = list())
# debug = TRUE

elbo0 <- mog_elbo(mog_normal_wishart, dgp$X, init$pars, init$fixed_pars, init$iNDK,
  amortization = list(features = NULL),
  debug = TRUE)

elbo0$elbo
#[1] -1031.255

# CAVI algorithm.

# family = mog_normal_wishart
# X = dgp$X
# K = dgp$K
# Xlags = NULL
# init_args = list(seed = 123321)
# maxiter = 500
# tol = 1e-8
# verbose = TRUE
# debug = FALSE

fit1 <- mog_cavi(mog_normal_wishart, dgp$X, dgp$K,
  init_args = list(seed = 123321),
  maxiter = 500, tol = 1e-8,
  verbose = TRUE, debug = FALSE)

"
iter =    1 | ELBO = -964.585131 | max mu = 2.0598
iter =    2 | ELBO = -944.086308 | max mu = 2.0662
iter =    3 | ELBO = -937.367517 | max mu = 2.0444
iter =    4 | ELBO = -935.006528 | max mu = 2.0287
iter =    5 | ELBO = -934.450343 | max mu = 2.0205
iter =    6 | ELBO = -934.346453 | max mu = 2.0167
iter =    7 | ELBO = -934.328261 | max mu = 2.0151
iter =    8 | ELBO = -934.325122 | max mu = 2.0144
iter =    9 | ELBO = -934.324583 | max mu = 2.0141
iter =   10 | ELBO = -934.324490 | max mu = 2.0140
iter =   11 | ELBO = -934.324475 | max mu = 2.0140
iter =   12 | ELBO = -934.324472 | max mu = 2.0140
"

##

init_aux <- init
init_aux$pars$mu_vs <- exp(init_aux$pars$mu_vs)
elbo0 <- mog2_elbo(x = x, pars = init_aux$pars,
  fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = forward_network, args = list(return_cache = FALSE)),
  debug = TRUE)

#elbo0$elbo
#[1] -1025.835

# Numerical minimization of the negative ELBO.

fit <- stats::optim(par = flatten_pars(init$pars), fn = mog2_neg_elbo_fn,
  model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = forward_network, args = list(return_cache = FALSE)),
  method = "BFGS", hessian = FALSE)

print(fit)
""

x_grid <- matrix(seq(min(x), max(x), length.out = 200), ncol = 1)

scores <- matrix(0, nrow = nrow(x_grid), ncol = dgp$K)

for (k in seq_len(dgp$K - 1)) {
  scores[, k + 1] <- forward_network(
    X = x_grid, etas = fitted_pars$etas[[k]],
    return_cache = FALSE)$Y
}

plot(x_grid[,1], scores[,2], type = "l")
abline(lm(scores[,2] ~ x_grid[,1]), col = 2)

plot(x_grid[,1], scores[,3], type = "l")
abline(lm(scores[,3] ~ x_grid[,1]), col = 2)
fitted_pars$mu_ms


fitted_pars
""

# Assess whether the neural-network learned a linear form.

Maux <- matrix(0.0, nrow = dgp$n, ncol = dgp$K)
for (k in seq_len(dgp$K-1))
  Maux[,k+1] <- forward_network(X = am_feats, etas = fitted_pars$etas[[k]], return_cache = FALSE)$Y

gammas <- softmax_rows(Maux)
#rowSums(gammas)

x_grid <- matrix(seq(min(x), max(x), length.out = 200), ncol = 1)
scores <- matrix(0, nrow = nrow(x_grid), ncol = dgp$K)
for (k in seq_len(dgp$K - 1))
  scores[,k+1] <- forward_network(X = x_grid, etas = fitted_pars$etas[[k]], return_cache = FALSE)$Y

#png("nn-shape.png")
op <- par(mfrow=c(2,1), mar=c(2,2,2,2))
plot(x_grid[,1], scores[,2], type = "l")
abline(lm(scores[,2] ~ x_grid[,1]), col = 2)

plot(x_grid[,1], scores[,3], type = "l")
abline(lm(scores[,3] ~ x_grid[,1]), col = 2)
par(op)
#dev.off()
