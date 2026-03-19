library(Matrix)
library(MetricGraph)
library(viridis)

# get the functions from the Rmd files
knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R"))
knitr::purl(here::here("matern_functions.Rmd"), output = here::here("matern_functions.R"))

# load the functions
source(here::here("all_functions.R"))
source(here::here("matern_functions.R"))

# build a graph with a mesh
graph_initial <- gets.graph.tadpole()
graph <- graph_initial$clone()

graph_initial$build_mesh(h = 1)
graph_initial$compute_fem()

mesh_loc <- graph_initial$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(y = 1) %>% 
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(mesh_loc,
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE, 
                       clear_obs = TRUE)

# add observations as vertices
graph$observation_to_vertex()
# get edge lengths
EDGE_LENGTHS <- graph$edge_lengths

# parameters
kappa <- 10
sigma <- 1
alpha <- 1.3
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 100
c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
# get rational approximation coefficients
coeff <- rSPDE:::interp_rational_coefficients(
  order = m, 
  type_rational_approx = "chebfun", 
  type_interp = "spline", 
  alpha = alpha)

r <- coeff$r; p <- coeff$p; k <- coeff$k

Qtilde_i <- list() 

for(order in 0:m){
  Qtilde_i[[paste0("m=",order)]] <- list()
  if(order == 0){
    P <- order
    ALPHA <- floor(alpha)
    r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
    correction_term <- rbind(cbind(r00_inverse, matrix(0, floor(alpha), floor(alpha))),
                             cbind(matrix(0, floor(alpha), floor(alpha)), r00_inverse))
  } else {
    P <- p[order]
    ALPHA <- alpha
    r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
    correction_term <- rbind(cbind(r00_inverse, matrix(0, ceiling(alpha), ceiling(alpha))),
                             cbind(matrix(0, ceiling(alpha), ceiling(alpha)), r00_inverse))
  }
  for(e in 1:length(EDGE_LENGTHS)){
    l_e <- EDGE_LENGTHS[e]
    if(order == 0){
      LOC <- c(0,l_e)
      DIVIDER <- k*kappa^2
      FACTOR <- 2*kappa*tau^2
    } else {
      LOC <- c(l_e, 0)
      DIVIDER <- r[order] * kappa^4
      FACTOR <- 2*c_alpha*sqrt(pi)*tau^2*kappa^3
      }
    aux <- matern.p.precision(loc = LOC, 
                              kappa = kappa, 
                              p = P,
                              equally_spaced = FALSE, 
                              alpha = ALPHA)
    Qtilde_i[[paste0("m=",order)]][[e]] <- (aux$Q - 0.5 * correction_term)*FACTOR
  }
  Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])/DIVIDER
}

Q1 <- Qtilde_i[[paste0("m=",0)]]

cbmat <- conditioning(graph,alpha=1)
nc1 <- 1:length(cbmat$S)
W <- Diagonal(dim(Q1)[1])
W <- W[,-nc1]
Qtilde1_uu <- t(W) %*% t(cbmat$T) %*% Q1 %*% (cbmat$T) %*% W 
index.obs1 <- gives.indices(graph = graph, factor = 2, constant = 2)
A0 <- cbmat$T[index.obs1,-nc1]

Qtilde_i <- Qtilde_i[-1] # remove m=0

graph$buildC(alpha = 2)
nc2 <- 1:length(c(1,graph$CoB$S))
TT <- graph$CoB$T
Tcomplete <- t(TT)[, c(ncol(TT), 1:(ncol(TT)-1))]

W2<-Diagonal(2*ceiling(alpha)*graph$nE)
W2 <- W2[,-nc2]

Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) t(W2) %*% t(Tcomplete) %*% Q %*% Tcomplete %*% W2) 

index.obs2 <- gives.indices(graph = graph, factor = 4, constant = 3)
Ai <- Tcomplete[index.obs2, -nc2]

A <- cbind(A0, do.call(cbind, rep(list(Ai), m)))
Qtilde_UU <- bdiag(Qtilde1_uu, do.call(bdiag, Qtilde_i_star_UU))*kappa^(2*alpha)

Sigma <- A %*% solve(Qtilde_UU, t(A)) 

True_Sigma <- gets_true_cov_mat(graph = graph_initial, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)


q <- graph_initial$plot_function(True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
graph_initial$plot_function(Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)

L_2_error = sqrt(as.double(t(graph_initial$mesh$weights)%*%(True_Sigma - Sigma)^2%*%graph_initial$mesh$weights))
print(L_2_error)

# function that translates from Vaibhav's code

gets_cov_mat_rat_approx_alpha_1_to_2 <- function(graph, kappa, tau, alpha, m){
  
  # get rational approximation coefficients
  coeff <- rSPDE:::interp_rational_coefficients(
    order = m, 
    type_rational_approx = "chebfun", 
    type_interp = "spline", 
    alpha = alpha)
  
  r <- coeff$r
  p <- coeff$p
  k <- coeff$k
  
  # compute constant c_alpha
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  
  # get edge lengths
  L_e <- graph$edge_lengths
  
  Qtilde_i <- list() 
  for(order in 0:m){
    if(order == 0){
      P <- order
      ALPHA <- floor(alpha)
      FACTOR <- (2*tau^2)/(k*kappa)
      r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
      correction_term <- rbind(cbind(r00_inverse, matrix(0, floor(alpha), floor(alpha))),
                               cbind(matrix(0, floor(alpha), floor(alpha)), r00_inverse))
    } else {
      P <- p[order]
      ALPHA <- alpha
      FACTOR <- (2*c_alpha*sqrt(pi)*tau^2)/(r[order] * kappa)
      r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
      correction_term <- rbind(cbind(r00_inverse, matrix(0, ceiling(alpha), ceiling(alpha))),
                               cbind(matrix(0, ceiling(alpha), ceiling(alpha)), r00_inverse))
    }
    Qtilde_i[[paste0("m=",order)]] <- list()
    for(e in 1:length(L_e)){
      Q_e <- matern.p.precision(loc = c(0, L_e[e]), #if (order == 0) c(0, L_e[e]) else c(L_e[e], 0), 
                                kappa = kappa, 
                                p = P,
                                equally_spaced = TRUE, 
                                alpha = ALPHA)$Q
      Qtilde_i[[paste0("m=",order)]][[e]] <- (Q_e - 0.5 * correction_term)*FACTOR*kappa^(2*alpha)
    }
    Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])
  }
  
  Qtilde_0 <- Qtilde_i[[paste0("m=",0)]] # extract Qtilde_0
  print(dim(Qtilde_0))
  Qtilde_i <- Qtilde_i[-1] # remove Qtilde_0
  
  # graph$.__enclos_env__$private$A()
  #####################################
  ## CASE m = 0
  #####################################
  COND_0 <- conditioning(graph = graph, alpha = 1)
  index.obs_0 <- gives.indices(graph = graph, factor = 2, constant = 2)
  nc_0 <- 1:length(COND_0$S) # number of constraints
  T_0 <- COND_0$T # change of basis matrix
  W_0 <- Diagonal(2*floor(alpha)*graph$nE)[,-nc_0] # matrix to remove constraints
  Qtilde_0_star_UU <- t(W_0) %*% t(T_0) %*% Qtilde_0 %*% (T_0) %*% W_0 
  A0 <- T_0[index.obs_0, -nc_0] # observation matrix after conditioning
  
  #####################################
  ## CASE m > 0
  #####################################
  graph$buildC(alpha = 2)
  COND_i <- graph$CoB
  index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  nc_i <- 1:length(c(1,COND_i$S)) # number of constraints
  T_i <- COND_i$T # change of basis matrix
  T_i <- t(T_i)[, c(ncol(T_i), 1:(ncol(T_i) - 1))] # column reordering
  W_i <- Diagonal(2*ceiling(alpha)*graph$nE)[,-nc_i] # matrix to remove constraints
  Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) t(W_i) %*% t(T_i) %*% Q %*% T_i %*% W_i)
  Ai <- T_i[index.obs_i, -nc_i] # observation matrix after conditioning
  
  
  #####################################
  ## Build matrix A and Q_UU
  #####################################
  A <- cbind(A0, do.call(cbind, rep(list(Ai), m)))
  Q_UU <- bdiag(Qtilde_0_star_UU, do.call(bdiag, Qtilde_i_star_UU))
  # Return Sigma
  Sigma <- A %*% solve(Q_UU, t(A)) 
  return(Sigma)
}

# this version is already ok, the next I just added with purrr
gets_cov_mat_rat_approx_alpha_1_to_2 <- function(graph, kappa, tau, alpha, m){
  
  # get rational approximation coefficients
  coeff <- rSPDE:::interp_rational_coefficients(
    order = m, 
    type_rational_approx = "chebfun", 
    type_interp = "spline", 
    alpha = alpha)
  
  r <- coeff$r
  p <- coeff$p
  k <- coeff$k
  
  # reparameterization
  nu <- alpha - 1/2
  sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  c_1 <- gamma(floor(alpha))/gamma(floor(alpha) - 0.5)
  # get edge lengths
  L_e <- graph$edge_lengths
  
  Qtilde_i <- list() 
  for(order in 1:m){
    r00_inverse <- solve(matern.p.joint(s = 0, 
                                        t = 0, 
                                        kappa = kappa, 
                                        p = p[order], 
                                        alpha = alpha))
    
    correction_term <- rbind(cbind(r00_inverse, matrix(0, ceiling(alpha), ceiling(alpha))),
                             cbind(matrix(0, ceiling(alpha), ceiling(alpha)), r00_inverse))
    Qtilde_i[[paste0("m=",order)]] <- list()
    for(e in 1:length(L_e)){
      Q_e <- matern.p.precision(loc = c(0, L_e[e]),
                                kappa = kappa, 
                                p = p[order],
                                equally_spaced = TRUE, 
                                alpha = alpha)$Q
      Qtilde_i[[paste0("m=",order)]][[e]] <- Q_e - 0.5 * correction_term
    }
    Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])*(2*kappa^(2*alpha)*c_alpha*sqrt(pi)*tau^2)/(r[order] * kappa)
  }
  #####################################
  ## CASE i = 0
  #####################################
  Qtilde_0_star_UU <- MetricGraph:::Qalpha1(theta = c(tau, kappa), 
                                            graph = graph, 
                                            BC = 3000, 
                                            build = TRUE)*c_1/(2 * k * c_alpha * kappa * sigma^2 * tau^2)
  A0 <- graph$.__enclos_env__$private$A()
  #####################################
  ## CASE i = 1,...,m
  #####################################
  graph$buildC(alpha = 2, edge_constraint = TRUE)
  COND_i <- graph$CoB
  index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  n_const <- length(COND_i$S)
  ind.const <- c(1:n_const)
  Tc <- COND_i$T[-ind.const, ]
  Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) Tc %*% Q %*% t(Tc)) 
  Ai <- t(Tc)[index.obs_i, ] # observation matrix after conditioning
  
  #####################################
  ## Build matrix A and Q_UU
  #####################################
  A <- cbind(A0, do.call(cbind, rep(list(Ai), m)))
  Q_UU <- bdiag(Qtilde_0_star_UU, do.call(bdiag, Qtilde_i_star_UU))
  # Return Sigma
  Sigma <- A %*% solve(Q_UU, t(A)) 
  return(Sigma)
}
  