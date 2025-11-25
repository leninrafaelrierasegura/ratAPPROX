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
sigma <- 1.4
alpha <- 1.8
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 100
Acc <- gamma(nu) / (kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2))
sigma2 <- Acc / tau^2
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
    r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = 0, alpha = floor(alpha)))
    correction_term <- rbind(cbind(r00_inverse, matrix(0, floor(alpha), floor(alpha))),
                             cbind(matrix(0, floor(alpha), floor(alpha)), r00_inverse))
  } else {
    r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = kappa^2*p[order], alpha = alpha))
    correction_term <- rbind(cbind(r00_inverse, matrix(0, ceiling(alpha), ceiling(alpha))),
                             cbind(matrix(0, ceiling(alpha), ceiling(alpha)), r00_inverse))
  }
  print(r00_inverse)
  for(e in 1:length(EDGE_LENGTHS)){
    l_e <- EDGE_LENGTHS[e]
    aux <- matern.p.precision(loc = c(0,l_e), 
                              kappa = kappa, 
                              p = ifelse(order == 0, order, p[order]),
                              equally_spaced = FALSE, 
                              alpha = ifelse(order == 0, floor(alpha), alpha))
    Qtilde_i[[paste0("m=",order)]][[e]] <- (aux$Q - 0.5 * correction_term) #/ ifelse(order == 0, k * sigma^2, r[order] * sigma^2)
    
  }
  Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])
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
nc2 <- 1:length(graph$CoB$S)
Tcomplete <- graph$CoB$T 

W2<-Diagonal(2*ceiling(alpha)*graph$nE)
W2 <- W2[,-nc2]

Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) t(W2) %*% t(Tcomplete) %*% Q %*% Tcomplete %*% W2) 

index.obs2 <- gives.indices(graph = graph, factor = 4, constant = 3)
Ai <- Tcomplete[index.obs2, -nc2]

A <- cbind(A0, do.call(cbind, rep(list(Ai), m)))
Qtilde_UU <- bdiag(Qtilde1_uu, do.call(bdiag, Qtilde_i_star_UU))

Sigma <- A %*% solve(Qtilde_UU, t(A)) 

True_Sigma <- gets_true_cov_mat(graph = graph_initial, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)


p <- graph_initial$plot_function(True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
graph_initial$plot_function(Sigma[,2], p = p, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)

L_2_error = sqrt(as.double(t(graph_initial$mesh$weights)%*%(True_Sigma - Sigma)^2%*%graph_initial$mesh$weights))
print(L_2_error)


matern.p.deriv(s=0,t=0,kappa=kappa,p=-34,alpha=2,deriv = 0)*sigma2
  
  