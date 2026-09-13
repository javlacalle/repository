## Bayesian variational inference by means of
## coordinate-ascent variational inference algorithm (CAVI)
## in the mixture of Gaussians (MoG) model described in the reference paper.
##
## Reference paper:
## Blei, D. M., Kucukelbir, A., and McAuliffe, J. D. (2017).
## "Variational Inference: A Review for Statisticians."
## In "Journal of the American Statistical Association", vol. 112, no. 518, pp. 859–877.
## https://doi.org/10.1080/01621459.2017.1285773


library(AAIMF)

dgp <- mog_sim(n = 500, K = 3, mu_sigma2 = 1, seed = 123)

print(dgp)

#head(dgp$x)
#head(dgp$c)
#head(dgp$C)

x <- dgp$x
#plot(x)

#init <- mog_init(x = x, m = 2, seed = 123321)
init <- mog_init(x = x, K = 2,
  mu_sigma2 = 1,
  amortization_features = NULL,
  seed = 123321)

init <- mog_init(x = x, K = 2,
  mu_sigma2 = 1,
  amortization_features = cbind(1, x), eta_sd = 0.01,
  seed = 123321)

model_struct <- get_pars_structure(init$pars)
vec_pars_to_list(flatten_pars(init$pars), model_struct)
flatten_pars(init$pars)

elbo0 <- mog_elbo_v0(x = x, pars = init$pars, fixed_pars = init$fixed_pars, debug = TRUE)
#elbo0

elbo <- mog_elbo(x = x, pars = init$pars, fixed_pars = init$fixed_pars, debug = TRUE)
#elbo

-mog_neg_elbo_fn(vec_pars = unpack(init$pars), model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars, amortization_features = NULL)

stats::optim(par = unpack(init$pars), fn = mog_neg_elbo_fn,
  model_struct = model_struct,
  x = x, fixed_pars = init$fixed_pars, amortization_features = NULL,
  method = "CG", hessian = FALSE) # method = "BFGS", "CG"


fit <- mog_cavi(x = x, m = 3,
    seed_init = 123321,
    maxiter = 500, tol = 1e-8,
    verbose = TRUE)

print(fit)

"""
> print(dgp)
Simulation from the Blei--Kucukelbir--McAuliffe Gaussian mixture model
n: 500
K: 3
sigma2: 1

Component means:
[1] -0.5604756 -0.2301775  1.5587083

Class counts:

  1   2   3
160 174 166

> fit <- mog_cavi(x = x, m = 3, seed_init = 123321, maxiter = 500, tol = 1e-8, verbose = TRUE)

iter =    1 | ELBO = -877.345737 | max mu = 1.6333
iter =    2 | ELBO = -862.007364 | max mu = 1.4830
iter =    3 | ELBO = -859.962592 | max mu = 1.4274
iter =    4 | ELBO = -859.616715 | max mu = 1.4069
iter =    5 | ELBO = -859.530903 | max mu = 1.4004
iter =    6 | ELBO = -859.494659 | max mu = 1.3993
iter =    7 | ELBO = -859.472239 | max mu = 1.4004
iter =    8 | ELBO = -859.456216 | max mu = 1.4021
iter =    9 | ELBO = -859.444280 | max mu = 1.4039
iter =   10 | ELBO = -859.435273 | max mu = 1.4056
iter =   11 | ELBO = -859.428433 | max mu = 1.4070
iter =   12 | ELBO = -859.423217 | max mu = 1.4083
iter =   13 | ELBO = -859.419222 | max mu = 1.4095
iter =   14 | ELBO = -859.416154 | max mu = 1.4104
iter =   15 | ELBO = -859.413788 | max mu = 1.4113
iter =   16 | ELBO = -859.411961 | max mu = 1.4120
iter =   17 | ELBO = -859.410544 | max mu = 1.4126
iter =   18 | ELBO = -859.409445 | max mu = 1.4132
iter =   19 | ELBO = -859.408589 | max mu = 1.4137
iter =   20 | ELBO = -859.407922 | max mu = 1.4141
iter =   21 | ELBO = -859.407401 | max mu = 1.4145
iter =   22 | ELBO = -859.406994 | max mu = 1.4148
iter =   23 | ELBO = -859.406675 | max mu = 1.4151
iter =   24 | ELBO = -859.406425 | max mu = 1.4153
iter =   25 | ELBO = -859.406229 | max mu = 1.4155
iter =   26 | ELBO = -859.406076 | max mu = 1.4157
iter =   27 | ELBO = -859.405955 | max mu = 1.4159
iter =   28 | ELBO = -859.405860 | max mu = 1.4160
iter =   29 | ELBO = -859.405785 | max mu = 1.4162
iter =   30 | ELBO = -859.405726 | max mu = 1.4163
iter =   31 | ELBO = -859.405679 | max mu = 1.4164
iter =   32 | ELBO = -859.405643 | max mu = 1.4165
iter =   33 | ELBO = -859.405614 | max mu = 1.4166
iter =   34 | ELBO = -859.405591 | max mu = 1.4167
iter =   35 | ELBO = -859.405573 | max mu = 1.4167
iter =   36 | ELBO = -859.405559 | max mu = 1.4168
iter =   37 | ELBO = -859.405548 | max mu = 1.4168
iter =   38 | ELBO = -859.405539 | max mu = 1.4169
iter =   39 | ELBO = -859.405532 | max mu = 1.4169

> print(fit)

Variational Bayes fit for mixture-of-Gaussians model
----------------------------------------------------------------
Converged: TRUE
Iterations: 39
Final ELBO: -859.4055

Model dimensions
----------------------------------------------------------------
Number of observations: 500
Number of mixture components m: 3

q(mu_s): Gaussian variational parameters
----------------------------------------------------------------
 component   mean_m_s variance_s2_k         sd
         1 -0.6651988   0.005963767 0.07722543
         2  0.1028175   0.005960177 0.07720218
         3  1.4169103   0.005968707 0.07725741

q(S): responsibility summary
----------------------------------------------------------------
 component expected_count mean_responsibility max_responsibility
         1       166.6793           0.3333585          0.8921378
         2       166.7803           0.3335605          0.4616523
         3       166.5405           0.3330810          0.9858386

"""
