## Variational Bayes / CAVI for the Generalized Autoregressive model
## with mixture-of-Gaussians innovations.
##
## Reference:
## Roberts and Penny (2002),
## "Variational Bayes for Generalized Autoregressive Models",
## IEEE Transactions on Signal Processing.
## https://doi.org/10.1109/TSP.2002.801921

library(AAIMF)

dgp <- ar_mog_sim(n = 500,
  model = list(ar = c(0.55, -0.25), means = c(0, 3), precisions = c(4, 1), probs = c(0.75, 0.25)),
  n.start = NA,
  seed = 123)

x <- dgp$x
#dgp$states
#plot(x)

init <- ar_mog_init(x = x, m = 2, p = 2, seed = 123)

elbo <- ar_mog_elbo(x, init)
#elbo

fit <- ar_mog_cavi(x = x, m = 2, p = 2,
    seed_init = 123,
    maxiter = 500, tol = 1e-8,
    verbose = TRUE)

print(fit)

"""
iter =    1 | ELBO = -770.885575 | max mu = 2.8911 | Eq_alpha = 0.8157
iter =    2 | ELBO = -729.031726 | max mu = 2.9358 | Eq_alpha = 5.1176
iter =    3 | ELBO = -723.866944 | max mu = 2.9813 | Eq_alpha = 5.3111
iter =    4 | ELBO = -722.631794 | max mu = 3.0223 | Eq_alpha = 5.4071
iter =    5 | ELBO = -722.095307 | max mu = 3.0579 | Eq_alpha = 5.4624
iter =    6 | ELBO = -721.761019 | max mu = 3.0874 | Eq_alpha = 5.4962
iter =    7 | ELBO = -721.536852 | max mu = 3.1111 | Eq_alpha = 5.5181
iter =    8 | ELBO = -721.390471 | max mu = 3.1297 | Eq_alpha = 5.5335
iter =    9 | ELBO = -721.298606 | max mu = 3.1440 | Eq_alpha = 5.5449
iter =   10 | ELBO = -721.242814 | max mu = 3.1550 | Eq_alpha = 5.5535
iter =   11 | ELBO = -721.209744 | max mu = 3.1633 | Eq_alpha = 5.5601
iter =   12 | ELBO = -721.190483 | max mu = 3.1696 | Eq_alpha = 5.5651
iter =   13 | ELBO = -721.179408 | max mu = 3.1743 | Eq_alpha = 5.5689
iter =   14 | ELBO = -721.173098 | max mu = 3.1779 | Eq_alpha = 5.5717
iter =   15 | ELBO = -721.169528 | max mu = 3.1805 | Eq_alpha = 5.5738
iter =   16 | ELBO = -721.167519 | max mu = 3.1825 | Eq_alpha = 5.5755
iter =   17 | ELBO = -721.166391 | max mu = 3.1840 | Eq_alpha = 5.5767
iter =   18 | ELBO = -721.165761 | max mu = 3.1851 | Eq_alpha = 5.5775
iter =   19 | ELBO = -721.165409 | max mu = 3.1860 | Eq_alpha = 5.5782
iter =   20 | ELBO = -721.165213 | max mu = 3.1866 | Eq_alpha = 5.5787
iter =   21 | ELBO = -721.165104 | max mu = 3.1870 | Eq_alpha = 5.5791
iter =   22 | ELBO = -721.165044 | max mu = 3.1874 | Eq_alpha = 5.5794
iter =   23 | ELBO = -721.165010 | max mu = 3.1876 | Eq_alpha = 5.5796
iter =   24 | ELBO = -721.164991 | max mu = 3.1878 | Eq_alpha = 5.5797
iter =   25 | ELBO = -721.164981 | max mu = 3.1880 | Eq_alpha = 5.5798
iter =   26 | ELBO = -721.164975 | max mu = 3.1881 | Eq_alpha = 5.5799

Variational Bayes fit for AR mixture-of-Gaussians-noise model
----------------------------------------------------------------
Converged: TRUE
Iterations: 26
Final ELBO: -721.165

Model dimensions
----------------------------------------------------------------
Number of observations used: 498
AR order p: 2
Number of mixture components m: 2

q(pi): Dirichlet variational parameters
----------------------------------------------------------------
 component   lambda      E_pi
         1 391.2711 0.7825421
         2 108.7289 0.2174579

q(mu_s): Gaussian variational parameters
----------------------------------------------------------------
 component   mean_m_s precision_v_s         sd
         1 0.03284127     1456.3224 0.02620422
         2 3.18807106      123.4083 0.09001767

q(beta_s): Gamma(shape = a_s, rate = r_s) variational parameters
----------------------------------------------------------------
 component shape_a_s rate_r_s   E_beta E_log_beta
         1 195.14554 52.29953 3.731305  1.3141938
         2  53.87446 47.02478 1.145661  0.1266723

q(w): multivariate Gaussian variational parameters
----------------------------------------------------------------
 coefficient       mean         sd
          w1  0.5319594 0.01615770
          w2 -0.2418708 0.01626351

Covariance matrix Sigma_w
       lag1   lag2
lag1  3e-04 -2e-04
lag2 -2e-04  3e-04

q(alpha): Gamma(shape = a_alpha, rate = r_alpha) variational parameters
----------------------------------------------------------------
 shape_a_alpha rate_r_alpha E_alpha E_log_alpha
          1.01    0.1810039 5.57999    1.148351

q(S): responsibility summary
----------------------------------------------------------------
 component expected_count mean_responsibility max_responsibility
         1       390.2711           0.7836769          0.9999592
         2       107.7289           0.2163231          1.0000000
"""

#

ar_mog_objective_map_em(x, init, debug = TRUE)
