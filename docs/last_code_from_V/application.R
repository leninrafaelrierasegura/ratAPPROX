library(MetricGraph)
library(rSPDE)
library(sf)
library(xtable)
library(scales)
library(gridExtra)
library(ggmap)
library(qs)

source("auxiliary_functions.R")
source("likelihood.R")
source("crossvalidation.R")


graph_pems <-  metric_graph$new(edges = pems_repl$edges, longlat = TRUE)

graph_pems$add_observations(data = pems_repl$data, group = "repl", normalized=TRUE)

graph_pems$add_observations(data = graph_pems$mutate(y = y-mean(y)), clear_obs=TRUE)

fit_fem = qread("data_files/fiiting_fem_with_L-BFGS-B.qs")

sigma_e_start <- exp(fit_fem$start_values[1])
range_start <- exp(fit_fem$start_values[3])
nu_start <- 1
sigma_start <- 1
theta0 <- c(log(c(sigma_e_start, sigma_start, range_start,nu_start)))

m = 4

# Crossvalidation for the proposed method (Rational approximation)

est_model_rational_default = fit_rational(y ~ -1, theta0, graph=graph_pems, m = m, method = "L-BFGS-B")

sigma_e_est = exp(est_model_rational_default$par$par[1])
sigma_est = exp(est_model_rational_default$par$par[2])
range_est <- exp(est_model_rational_default$par$par[3])
# nu_est    = 1.5*exp(est_model_rational$par$par[4])/(1+exp(est_model_rational$par$par[4]))
 nu_est = 2.5*exp(est_model_rational_default$par$par[4])/(1+exp(est_model_rational_default$par$par[4]))
kappa_est = sqrt(8*nu_est)/range_est
tau_est <- sqrt(gamma(nu_est) / (sigma_est^2 * kappa_est^(2 * nu_est) * (4 * pi)^(1 / 2) * gamma(nu_est + 1 / 2)))

c(sigma_e_est, sigma_est, range_est, nu_est)

crossvalid_rational = crossvalidation_rational_covariates(est_model = est_model_rational_default, 
                                                          graph = graph_pems, m=m, factor =1)[["scores"]]


# Crossvalidation for Whittle-Matérn with alpha=1 and alpha=2

fit_alpha1 <- graph_lme(y ~ -1, graph=graph_pems, BC=0,        # In case of boundary correction change BC=1
            model = list(type = "WhittleMatern", alpha = 1))

fit_alpha2 <- graph_lme(y ~ -1, graph=graph_pems, BC=0,         # In case of boundary correction change BC=1
            model = list(type = "WhittleMatern", alpha = 2))

crossvalid_alpha1 = posterior_crossvalidation(fit_alpha1, factor=1)[["scores"]]
crossvalid_alpha2 = posterior_crossvalidation(fit_alpha2, factor=1)[["scores"]]


# Crossvalidation for covariance based rational SPDE approach (FEM)

fit_fem <- graph_lme(y ~ -1, model = "wm", graph=graph_pems)
crossvalid_fem = posterior_crossvalidation(fit_fem, factor=1)[["scores"]]


# Crossvalidation for isotropic exponential covariance model by Anderes et al.

fit_isoexp <- graph_lme(y ~ -1, graph=graph_pems, 
                model = list(type = "isoCov"))

crossvalid_isoexp = posterior_crossvalidation(fit_isoexp, factor=1)[["scores"]]