library(rSPDE)
library(rdist)
library(MetricGraph)
library(Matrix)
library(sp)
library(tidyr)

source("auxiliary_functions.R")

CRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) + 0.5 * Exx(mu, sigma))
}

#' @noRd
SCRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) / Exx(mu, sigma) - 0.5 * log(Exx(mu, sigma)))
}


#' @noRd
LS <- function(y, mu, sigma)
{
  return(dnorm(y, mean = mu, sd = sigma, log = TRUE))
}


#' @noRd
Exx <- function(mu, sigma) {
  #X-X' = N(0,2*sigma^2)
  return(Efnorm(0, sqrt(2) * sigma))
}



#compute E[|X-y|] when X is N(mu,sigma^2)

#' @noRd
Exy <- function(mu, sigma, y) {
  #X-y = N(mu-y,sigma^2)
  return(Efnorm(mu - y, sigma))
}

#' @noRd
Efnorm <- function(mu, sigma) {
  return(sigma * sqrt(2 / pi) * exp(-(mu ^ 2) / (2 * sigma ^ 2)) + mu * (1 -
                                                                           2 * pnorm(-mu / sigma)))
}

crossvalidation_rational = function(graph, m, factor=1, kappa_est, tau_est, nu_est, sigma_e_est){

 y_graph = graph$.__enclos_env__$private$data[[as.character("y")]]

# data <- graph$.__enclos_env__$private$data

# y_term <- stats::terms(formula)[[2]]

# y_graph <- eval(y_term, envir = data, enclos = parent.frame())
# y_graph <- as.numeric(y_graph)

ind <- 1:length(y_graph)
repl_vec <- graph$.__enclos_env__$private$data[[".group"]]
repl <- unique(repl_vec)


graph$observation_to_vertex(mesh_warning = FALSE)


if (nu_est>0 && nu_est<0.5){
  Q_tmp <- precision_alpha1(kappa = kappa_est, tau = tau_est,#1/reciprocal_tau,
                      nu = nu_est, graph = graph, m = m)
} else if(nu_est>0.5 && nu_est<1.5){
Q_tmp <- precision_alpha2(kappa = kappa_est, tau = tau_est,#1/reciprocal_tau,
                      nu = nu_est, graph = graph, m = m)
} else {
  stop("not implemented")
}

Q = Q_tmp$Q
A = Q_tmp$A

Sigma = A %*% solve(Q, t(A))

Sigma.o <- Sigma
diag(Sigma.o) <- diag(Sigma.o) + sigma_e_est^2

# ord = match_order(graph$V, graph$coordinates(graph$get_PtE()))
# y_tmp = apply(y, 2, function(col) col[ord])
# y_tmp = matrix(as.vector(y_tmp), ncol=1)

# y_graph = y_tmp[,1]

n_obs <- sum(repl_vec == repl[1])

mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, n_obs)
mae <- rmse <- rep(0, n_obs)

 for(i in 1:n_obs){

  idx_repl <- repl_vec == repl[1]
  y_graph_repl <- y_graph[idx_repl]

  y_cv <- y_graph_repl[-i]
  v_cv <- y_cv
  
   mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], v_cv)
   Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                              Sigma.o[-i, i])
   var.p[i] <- diag(Sigma.p)
  
  logscore[i] <- LS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  crps[i] <- CRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  scrps[i] <- SCRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  mae[i] <- abs(y_graph_repl[i] - mu.p[i])
  rmse[i] <- (y_graph_repl[i] - mu.p[i])^2

   if(length(repl)>1){
        for (j in 2:length(repl)) {
          y_graph_repl <- y_graph[repl_vec == repl[j]]
          y_cv <- y_graph_repl[-i]
          v_cv <- y_cv

    mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], v_cv)
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                                Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
   
   logscore[i] <- logscore[i] + LS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   crps[i] <- crps[i] + CRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   scrps[i] <- scrps[i] + SCRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   mae[i] <- mae[i] + abs(y_graph_repl[i] - mu.p[i])
   rmse[i] <- rmse[i] + (y_graph_repl[i] - mu.p[i])^2
   }
  }
 }

res <- list(mu = mu.p,
            var = var.p,
        logscore = -factor * mean(logscore/length(repl), na.rm = TRUE),
        crps = -factor * mean(crps/length(repl), na.rm = TRUE),
        scrps = -factor * mean(scrps/length(repl), na.rm = TRUE),
        mae = factor * mean(mae/length(repl), na.rm = TRUE),
        rmse = factor * sqrt(mean(rmse/length(repl), na.rm = TRUE)))

 return(res)

}

crossvalidation_rational_covariates = function(est_model, graph, m, factor=1){

graph$observation_to_vertex(mesh_warning = FALSE)

sigma_e_est = exp(est_model$par$par[1])
sigma_est = exp(est_model$par$par[2])
range_est <- exp(est_model$par$par[3])
nu_est = 2.5*exp(est_model$par$par[4])/(1+exp(est_model$par$par[4]))

#  nu_est = exp(est_model$par$par[4])

kappa_est = sqrt(8*nu_est)/range_est
tau_est <- sqrt(gamma(nu_est) / (sigma_est^2 * kappa_est^(2 * nu_est) * (4 * pi)^(1 / 2) * gamma(nu_est + 1 / 2)))

if (nu_est>0 && nu_est<0.5){
  Q_tmp <- precision_alpha1(kappa = kappa_est, tau = tau_est,#1/reciprocal_tau,
                      nu = nu_est, graph = graph, m = m)
} else if(nu_est>0.5 && nu_est<1.5){
Q_tmp <- precision_alpha2(kappa = kappa_est, tau = tau_est,#1/reciprocal_tau,
                      nu = nu_est, graph = graph, m = m)
} else if (nu_est>1.5 && nu_est<1.5){
Q_tmp <- shiftedrational_Qalpha3(kappa = kappa_est, tau = tau_est,#1/reciprocal_tau,
                      nu = nu_est, graph = graph, m = m)

}else{
  stop("not implemented")
}

Q = Q_tmp$Q
A = Q_tmp$A

Sigma = A %*% solve(Q, t(A))

Sigma.o <- Sigma
diag(Sigma.o) <- diag(Sigma.o) + sigma_e_est^2


# data <- graph$.__enclos_env__$private$data

# y_term <- stats::terms(formula)[[2]]

# y_graph <- eval(y_term, envir = data, enclos = parent.frame())
# y_graph <- as.numeric(y_graph)

if(!is.matrix(est_model$model_matrix)){
    est_model$model_matrix <- matrix(est_model$model_matrix, ncol=1)
  }

y_graph = est_model$model_matrix[,1]
ind <- 1:length(y_graph)

 if(ncol(est_model$model_matrix) > 1){
    X_cov <- est_model$model_matrix[,-1]
  } else{
    X_cov <- NULL
  }

if (!is.null (X_cov)){
beta_cov = est_model$par$par[5:(4+ncol(X_cov))]
}


repl_vec <- graph$.__enclos_env__$private$data[[".group"]]
repl <- unique(repl_vec)

    # ord = match_order(graph$V, graph$coordinates(graph$get_PtE()))
    # y_tmp = apply(y, 2, function(col) col[ord])
    # y_tmp = matrix(as.vector(y_tmp), ncol=1)

    # y_graph = y_tmp[,1]

n_obs <- sum(repl_vec == repl[1])

mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, n_obs)
mae <- rmse <- rep(0, n_obs)

 for(i in 1:n_obs){

  idx_repl <- repl_vec == repl[1]
  y_graph_repl <- y_graph[idx_repl]

  y_cv <- y_graph_repl[-i]
  v_cv <- y_cv

  if(!is.null(X_cov)){
          X_cov_repl <- X_cov[idx_repl, , drop = FALSE]
          v_cv <- v_cv - as.vector(X_cov_repl[-i, , drop = FALSE] %*% beta_cov)
          mu_fe <- as.vector(X_cov_repl[i, , drop = FALSE] %*% beta_cov)
        } else {
          mu_fe <- 0
    }
  
   mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], v_cv) + mu_fe
   Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                              Sigma.o[-i, i])
   var.p[i] <- diag(Sigma.p)
  
  logscore[i] <- LS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  crps[i] <- CRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  scrps[i] <- SCRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
  mae[i] <- abs(y_graph_repl[i] - mu.p[i])
  rmse[i] <- (y_graph_repl[i] - mu.p[i])^2

   if(length(repl)>1){
        for (j in 2:length(repl)) {
          y_graph_repl <- y_graph[repl_vec == repl[j]]
          y_cv <- y_graph_repl[-i]
          v_cv <- y_cv

          if(!is.null(X_cov)){
            X_cov_repl <- X_cov[idx_repl,, drop = FALSE]
            v_cv <- v_cv - as.vector(X_cov_repl[-i, , drop = FALSE] %*% beta_cov)
            mu_fe <- as.vector(X_cov_repl[i, , drop = FALSE] %*% beta_cov)
          } else {
            mu_fe <- 0
          }

    mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], v_cv) + mu_fe
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                                Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
   
   logscore[i] <- logscore[i] + LS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   crps[i] <- crps[i] + CRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   scrps[i] <- scrps[i] + SCRPS(y_graph_repl[i], mu.p[i], sqrt(var.p[i]))
   mae[i] <- mae[i] + abs(y_graph_repl[i] - mu.p[i])
   rmse[i] <- rmse[i] + (y_graph_repl[i] - mu.p[i])^2
   }
  }
 }

res <- list(mu = mu.p,
            var = var.p,
        logscore = -factor * mean(logscore/length(repl), na.rm = TRUE),
        crps = -factor * mean(crps/length(repl), na.rm = TRUE),
        scrps = -factor * mean(scrps/length(repl), na.rm = TRUE),
        mae = factor * mean(mae/length(repl), na.rm = TRUE),
        rmse = factor * sqrt(mean(rmse/length(repl), na.rm = TRUE)))

 return(res)

}



# Test

# range <- 0.2
# sigma <- 1.3
# sigma_e <- 0.1
# nu = 0.28
# # kappa = 3
# kappa = sqrt(8*nu)/range
# m = 2

# n_repl = 20

# line1 <- Line(rbind(c(0,0),c(1,0)))
# line2 <- Line(rbind(c(0,1/(1+pi/4)),c(0,0)))
# line3 <- Line(rbind(c(-1/(1+pi/4),1/(1+pi/4)),c(0,1/(1+pi/4))))
# theta <- seq(from=pi,to=3*pi/2,length.out = 50)
# line4 <- Line(cbind(sin(theta)/(1+pi/4),(1+ cos(theta))/(1+pi/4)))
# Lines = SpatialLines(list(Lines(list(line1),ID="1"),
#                           Lines(list(line2),ID="2"),
#                           Lines(list(line3),ID="3"),
#                           Lines(list(line4),ID="4")))

# graph <- metric_graph$new(lines = Lines)

# graph$build_mesh(n = 100)

# # graph$plot()

# op_cov <- matern.operators(
#   graph = graph, nu = nu,
#   range = range, sigma = sigma, d = 1, m = m,
#   parameterization = "matern", type_rational_approximation = "brasil"
# )

# u = simulate(op_cov, nsim = n_repl)
# y <- u + sigma_e*matrix(rnorm(length(graph$mesh$VtE[,1]) * n_repl), ncol = n_repl)

# PtE <- graph$mesh$VtE # Get mesh locations

# df_graph <- data.frame(y = as.matrix(y), edge_number = PtE[,1],
#                         distance_on_edge = PtE[,2])


# df_graph <- pivot_longer(df_graph, cols = `y.1`:`y.20`, names_to = "repl", values_to = "y")
# graph$add_observations(data = df_graph, normalized = TRUE, group = "repl", clear_obs=TRUE)

# fit <- graph_lme(y~-1, model = "wm", graph=graph)

# # krig_fem = predict(fit, newdata = data.frame(edge_number = PtE[,1],
# #                                       distance_on_edge = PtE[,2]), normalized = TRUE, return_as_list = TRUE)

# crossvalid_fem = posterior_crossvalidation(fit, factor=1000)[["scores"]]


# graph2 = graph$clone()
# graph_tmp = graph2$clone()

# graph$observation_to_vertex()

# df_graph2 <- data.frame(edge_number = graph2$mesh$VtE[,1],         # Prediction locations
#                         distance_on_edge = graph2$mesh$VtE[,2])

# sigma_e_start <- exp(fit$start_values[1])
# range_start <- exp(fit$start_values[3])
# nu_start <- 1
# sigma_start <- 1

# theta0 <- c(log(c(sigma_e_start, sigma_start, range_start, nu_start)))

#     likelihood_new <- function(theta, graph, y, m, repl, BC, parameterization) {
#       l_tmp <- tryCatch(likelihood_rational(theta, graph=graph, y=y, m=m, repl=repl, BC=BC, parameterization = parameterization),
#         error = function(e) {
#           return(NULL)
#         }
#       )
#       if (is.null(l_tmp)) {
#         return(10^100)
#       }
#       return(l_tmp)
#     }

# par <- optim(theta0, likelihood_new, method = "L-BFGS-B", graph = graph, y = y, m = m, repl= NULL, BC=0, parameterization = "matern")

#  sigma_e_est = exp(par$par[1])
# # sigma_e_est = fit$coeff$measurement_error
# sigma_est = exp(par$par[2])
# range_est <- exp(par$par[3])
# nu_est    = 1.5*exp(par$par[4])/(1+exp(par$par[4]))
# kappa_est = sqrt(8*nu_est)/range_est
# tau_est <- sqrt(gamma(nu_est) / (sigma_est^2 * kappa_est^(2 * nu_est) * (4 * pi)^(1 / 2) * gamma(nu_est + 1 / 2)))


# crossvalid = crossvalidation_rational(graph = graph2, m, factor=1000, y, kappa_est, tau_est, nu_est, sigma_e_est)



# object = fit
# if(!is.matrix(object$model_matrix)){
#     object$model_matrix <- matrix(object$model_matrix, ncol=1)
#   }

#   yy <- object$model_matrix[,1]

