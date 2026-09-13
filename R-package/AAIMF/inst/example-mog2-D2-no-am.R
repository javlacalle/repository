library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(c(-2, -1), c(2, 1)),
  Sigmas = list(diag(c(0.4, 0.8)), matrix(c(0.8, 0.3, 0.3, 0.6), 2, 2)),
  seed = 123)
#X <- dgp$X
#print(dgp)

dput(dgp, file="test-dgp.txt")

# build model

model <- make_model('mog_normal_wishart', dgp$X, dgp$K,
  priors = list())

# Initialize the parameters of the model.

# X = dgp$X
# K = dgp$K
# fixed_pars = list(alpha0 = 1, beta0 = 1, m0 = NULL, W0 = NULL, nu0 = NULL)
# amortization = list(features = NULL, type = c("linear", "network"), args = list())
# kmeans.nstart = 10
# kmeans.prob = 0.95
# seed = 123321

model1 <- mog_init(model, seed = 123321)
  #kmeans.nstart = 10, kmeans.prob = 0.95

model2 <- mog_nw_init(model)

for (label in names(model$priors))
{
  if (!identical(model1$priors[[label]], model2$priors[[label]]))
    print(label)
}

# ELBO

model <- model1

elbo <- mog_elbo(model,
  amortization = list(features = NULL, fn = NULL, args = list()),
  debug = TRUE)
elbo$elbo
# [1] -1628.641

mog_nw_elbo(model)$elbo
# [1] -1628.641

# CAVI

fit <- mog_cavi(model,
  maxiter = 500, tol = 1e-8,
  verbose = TRUE, debug = FALSE)

"
iter =    1 | ELBO = -1561.951492
iter =    2 | ELBO = -1539.838506
iter =    3 | ELBO = -1531.614296
iter =    4 | ELBO = -1529.001212
iter =    5 | ELBO = -1528.518598
iter =    6 | ELBO = -1528.450377
iter =    7 | ELBO = -1528.441670
iter =    8 | ELBO = -1528.440597
iter =    9 | ELBO = -1528.440466
iter =   10 | ELBO = -1528.440450
iter =   11 | ELBO = -1528.440448
"

#


##

#init$pars

dput(as_bishop_gmm_init_elbo(init), file="test-init.txt")

model_struct <- get_pars_structure(init$pars)

# check for consistency of parameters:
# mog2_init() should return parameter values in agreement with CAVI update expressions.
tmp <- mog2_update(dgp$X, init$pars$rprobs, init$fixed_pars, init$dim['K'], init$dim['D'])
for (nm in names(tmp))
  print(identical(tmp[[nm]], init$pars[[nm]]))
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE

# Evaluate the ELBO at the initial parameter values.

# X = dgp$X
# pars = init$pars
# fixed_pars = init$fixed_pars
# dim = init$dim
# amortization = list(features = NULL, fn = NULL, args = list())
# debug = TRUE

elbo0 <- mog2_elbo(dgp$X, init$pars, init$fixed_pars, init$dim,
  amortization = list(features = NULL),
  debug = TRUE)

elbo0$elbo
#[1] -1628.641

# CAVI algorithm.

# X = dgp$X
# K = dgp$K
# init_args = list(seed = 123321)
# maxiter = 500
# tol = 1e-8
# verbose = TRUE
# debug = FALSE

fit1 <- mog2_cavi(dgp$X, dgp$K,
  init_args = list(seed = 123321),
  maxiter = 500, tol = 1e-8,
  verbose = TRUE, debug = FALSE)

"
iter =    1 | ELBO = -1561.951492 | max mu = 2.0666
iter =    2 | ELBO = -1539.838506 | max mu = 2.0694
iter =    3 | ELBO = -1531.614296 | max mu = 2.0450
iter =    4 | ELBO = -1529.001212 | max mu = 2.0300
iter =    5 | ELBO = -1528.518598 | max mu = 2.0231
iter =    6 | ELBO = -1528.450377 | max mu = 2.0203
iter =    7 | ELBO = -1528.441670 | max mu = 2.0193
iter =    8 | ELBO = -1528.440597 | max mu = 2.0190
iter =    9 | ELBO = -1528.440466 | max mu = 2.0188
iter =   10 | ELBO = -1528.440450 | max mu = 2.0188
iter =   11 | ELBO = -1528.440448 | max mu = 2.0188
"
