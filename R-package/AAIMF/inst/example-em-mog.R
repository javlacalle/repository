library(AAIMF)

# Generate sample data.

dgp <- mog_sim(N = 500,
  pis = c(0.35, 0.65),
  mus = rbind(-2, 2), Sigmas = list(matrix(0.4, 1, 1), matrix(0.8, 1, 1)),
  seed = 123)

#

# NOTE if model$D == 1, there is no difference between 'mog_full_cov' and 'mog_diagonal_cov', right ???

model <- make_model(family = 'mog_full_cov', #'mog_full_cov' #'mog_diagonal_cov'
  dgp$X, dgp$K,
  #p = 0L, Xlags = NULL,
  priors = list(),
  amortization = list(features = NULL, fn = NULL, args_fn = list(), finit = NULL, init_args = list()))

model <- em_init(model,
  kmeans_nstart = 10, kmeans_prob = 0.95,
  seed = NULL)

em_loglik(model)
# before: probably including constant
# [1] -972.8883
# last:
# [1] -513.419

em_upd_resp(model, TRUE)

#em_mog_fullcov_upd(model)
#em_mog_diagcov_upd(model)

fit <- em_run(model,
  maxiter = 500, tol = 1e-8,
  verbose = TRUE, debug = FALSE)

"
iter =    1 | logLik = -932.567546 | max mu = 2.0675
iter =    2 | logLik = -919.663605 | max mu = 2.0727
iter =    3 | logLik = -914.738207 | max mu = 2.0496
iter =    4 | logLik = -913.313732 | max mu = 2.0330
iter =    5 | logLik = -913.028313 | max mu = 2.0244
iter =    6 | logLik = -912.976163 | max mu = 2.0204
iter =    7 | logLik = -912.966679 | max mu = 2.0186
iter =    8 | logLik = -912.964951 | max mu = 2.0178
iter =    9 | logLik = -912.964635 | max mu = 2.0174
iter =   10 | logLik = -912.964577 | max mu = 2.0173
iter =   11 | logLik = -912.964567 | max mu = 2.0172
iter =   12 | logLik = -912.964565 | max mu = 2.0172
"

fit

"
Expectation-Maximization fit for mixture-of-Gaussians model
Diagonal covariance matrix
====================================================

Model dimensions
----------------------------------------------------------------
Number of observations: 500
Number of variables D: 1
Number of mixture components K: 2

Convergence
----------------------------------------------------------------
Converged: TRUE
Iterations: 12
Initial Log-Likelihood: -972.8883
Final Log-Likelihood: -912.9646

Fitted parameters summary
----------------------------------------------------------------

Mixture probabilities:
[1] 0.6563 0.3437

Component means, m_k:
        [,1]
[1,]  2.0172
[2,] -1.9660

Component precisions:
       [,1]
[1,] 1.2458
       [,1]
[1,] 2.3869

Component variances:
       [,1]
[1,] 0.8027
      [,1]
[1,] 0.419
"

# confidence intervals
# obtained upon the inverse of the numerical Hessian

ci <- confint(fit, level = 0.95)
ci
"
                  2.5 %  97.5 %
pis[1]           0.6128  0.6972
pis[2]           0.3028  0.3872
mus[1,1]         1.9164  2.1180
mus[2,1]        -2.0686 -1.8633
Lambdas[1][1,1]  1.0502  1.4777
Lambdas[2][1,1]  1.8723  3.0428
"
