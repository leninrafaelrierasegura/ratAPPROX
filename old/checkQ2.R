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
alpha <- 1.8
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
      DIVIDER <- k * kappa^2
      FACTOR <- sqrt(pi)*2*tau^2*kappa^3
    } else {
      LOC <- c(l_e, 0)
      DIVIDER <- r[order] * kappa^4
      FACTOR <- c_alpha*sqrt(pi)*2*tau^2*kappa^3
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
Q1

  