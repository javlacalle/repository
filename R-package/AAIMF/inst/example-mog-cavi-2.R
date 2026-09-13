
library(AAIMF)

dgp <- mog_sim(n = 500, K = 3, mu_sigma2 = 1, seed = 123)
x <- dgp$x

#print(dgp)
"
Simulation from the Blei--Kucukelbir--McAuliffe Gaussian mixture model
n: 500
K: 3
sigma2: 1

Component means:
[1] -0.5604756 -0.2301775  1.5587083

Class counts:

  1   2   3
160 174 166
"

am_feats <- cbind(1,x) # cbind(1,x) # NULL # cbind(x)

init <- mog_init(x = x, K = dgp$K,
  mu_sigma2 = 1,
  #amortization_features = amortization_features,
  amortization = list(features = am_feats, type = "linear", args = list(sd=0.01)),
  seed = 123321)

# FIXME am_feats should have dimension 1, or input_dim required adjustment???

init <- mog_init(x = x, K = dgp$K,
  mu_sigma2 = 1,
  amortization = list(features = am_feats, type = "network",
                      args = list(input_dim = 1, width = 20, output_dim = 1, seed = NULL)),
  seed = 123321)

#etas <- do.call("init_network", list(input_dim = 1, width = 20, output_dim = 1, seed = NULL))
#etas <- do.call("init_network", list())

init$pars
"
$gammas
NULL

$mu_ms
[1] 0.2712656 0.2882381 0.2938538

$mu_vs
[1] -5.115422 -5.126046 -5.124433

$etas
[,1]        [,2]
[1,] 0.007964754 0.005404655
[2,] 0.009673793 0.012876865
"

model_struct <- get_pars_structure(init$pars)

init_aux <- init
if (!is.null(am_feats)) {
  # this adjustment is required for comparison of mog_neg_elbo_fn() below;
  # when amortization_features != NULL the amortization approach is followed by default,
  # so mog_init() returns log(mu_vx) rather than the standard scale.
  init_aux$pars <- vec_pars_to_list(flatten_pars(init$pars), model_struct, list("mu_vs", exp))
}

elbo0 <- mog_elbo(x = x, pars = init_aux$pars,
  fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = function(X, etas) X %*% etas),
  debug = TRUE)
#[1] -908.2315

init_aux <- init
init_aux$pars$mu_vs <- exp(init_aux$pars$mu_vs)
elbo0 <- mog_elbo(x = x, pars = init_aux$pars,
  fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = forward_network, args = list(return_cache = FALSE)),
  debug = TRUE)

# Use simple softmax regression (suggested by CAVI update):
# varphi_eta(x_i) = softmax(eta_{k0} + eta_{k1} x_i)

amortization_function <- function(X, etas) X %*% etas

elbo0 <- mog_elbo(x = x, pars = init_aux$pars, fixed_pars = init$fixed_pars,
  amortization = list(features = amortization_features, fn = amortization_function),
  debug = TRUE)

elbo1 <- mog_neg_elbo_fn(vec_pars = flatten_pars(init$pars), model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars, amortization_features = amortization_features)

all.equal(elbo0$elbo, -elbo1)
#TRUE


init_aux <- init
init_aux$pars$mu_vs <- exp(init_aux$pars$mu_vs)
elbo0 <- mog_elbo(x = x, pars = init_aux$pars,
  fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = function(X, etas) X%*% etas),
  debug = TRUE)
#elbo0$elbo
#[1] -908.2315


# NOTE Doing numerical optimization is sensible only with non-nulll amortization_features;
# otherwise several parameters (gammas increase with length(x)).
#
# Yet, amortization with a neural-network involves also a large number of parameters
# (see stochastic optimization/approximation as in Ann-VI-linear-model.R).

fit <- stats::optim(par = flatten_pars(init$pars), fn = mog_neg_elbo_fn,
  model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars,
  amortization = list(features = am_feats, fn = function(X, etas) X%*% etas),
  method = "BFGS", hessian = FALSE)


fit <- stats::optim(par = flatten_pars(init$pars), fn = mog_neg_elbo_fn,
  model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars, amortization_features = amortization_features,
  method = "BFGS", hessian = FALSE)

fit
"
$par
 [1] -0.6644295  1.4172172  0.1015497 -5.1224325 -5.1222743 -5.1211981
 [7] -0.7835135  2.0816308  0.2155832  0.7659691

$value
[1] 859.4055

$counts
function gradient
      87       32

$convergence
[1] 0

$message
NULL
"

fitted_pars <- vec_pars_to_list(fit$par, model_struct, list("mu_vs", exp))
fitted_pars$mu_ms
#[1] -0.6644295  1.4172172  0.1015497

print(fit)
"
$par
[1] -0.6644295  1.4172172  0.1015497 -5.1224325 -5.1222743 -5.1211981
[7] -0.7835135  2.0816308  0.2155832  0.7659691

$value
[1] 859.4055

$counts
function gradient
87       32

$convergence
[1] 0

$message
NULL
"

#vec_pars_to_list(fit$par, model_struct)
vec_pars_to_list(fit$par, model_struct, list("mu_vs", exp))
"
$gammas
numeric(0)

$mu_ms
[1] -0.6644295  1.4172172  0.1015497

$mu_vs
[1] 0.005961504 0.005962447 0.005968868

$etas
[,1]      [,2]
[1,] -0.7835135 0.2155832
[2,]  2.0816308 0.7659691
"
