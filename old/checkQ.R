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

graph_initial$build_mesh(h = 0.1)

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

# get rational approximation coefficients
coeff <- rSPDE:::interp_rational_coefficients(
  order = m, 
  type_rational_approx = "chebfun", 
  type_interp = "spline", 
  alpha = alpha)

r <- coeff$r
p <- coeff$p
k <- coeff$k

Qtilde_i <- list() # 60 x 60 # 12 x 12
A_mat_list <- list() # 2 x 60 # 2 x 12
A_list <- list()
for(order in 0:m){
  Qtilde_i[[paste0("m=",order)]] <- list()
  A_mat_list[[paste0("m=",order)]] <- list()
  for(e in 1:length(EDGE_LENGTHS)){
    l_e <- EDGE_LENGTHS[e]
    aux <- matern.p.precision(loc = c(0,l_e), 
                              kappa = kappa, 
                              p = ifelse(order == 0, order, p[order]),
                              equally_spaced = FALSE, 
                              alpha = alpha)
    Qtilde_i[[paste0("m=",order)]][[e]] <- aux$Q * ifelse(order == 0, k, r[order]) * sigma^2
    A_mat_list[[paste0("m=",order)]][[e]] <- aux$A
  }
  Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])
  A_mat_list[[paste0("m=",order)]] <- do.call(cbind, A_mat_list[[paste0("m=",order)]])
  A_list[[paste0("m=",order)]] <- aux$A
}

A <- do.call(cbind, A_list) # 2 x 12 # 2 x 12

graph$buildC(alpha = 2)
n_const <- length(graph$CoB$S)
ind.const <- c(1:n_const)
Tcomplete <- graph$CoB$T # 60 x 60 # 12 x 12
TUU <- Tcomplete[-ind.const, ]  # 31 x 60 # 7 x 12
Qtilde_i_star <- lapply(Qtilde_i, function(Q) Tcomplete %*% Q %*% t(Tcomplete)) # 60 x 60 # 12 x 12
Qtilde_i_star_UU <- lapply(Qtilde_i_star, function(Q) Q[-ind.const, -ind.const]) # 31 x 31 # 7 x 7
TUU_t_Qtilde_i_star_UU_TUU <- lapply(Qtilde_i_star_UU, function(Q) t(TUU) %*% Q %*% TUU) # 60 x 60 # 12 x 12
Q_UU <- bdiag(TUU_t_Qtilde_i_star_UU_TUU) # 180 x 180 # 36 x 36





A_tmp <- t(TUU)
index.obs1 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
idx_temp <- which(idx_temp)
return(idx_temp[1])})

index.obs1 <- (index.obs1-1)*4+1
index.obs2 <- NULL
na_obs1 <- is.na(index.obs1)
if(any(na_obs1)){
  idx_na <- which(na_obs1)
  PtV_NA <- graph$PtV[idx_na]
  index.obs2 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
  idx_temp <- which(idx_temp)
  return(idx_temp[1])})
  index.obs1[na_obs1] <- (index.obs2-1)*4 + 3                                                                      
}
A_tmp <- A_tmp[index.obs1,] # 15 x 31 # 3 x 7

aux_list <- lapply(Qtilde_i_star_UU, function(Q) A_tmp %*% solve(Q, t(A_tmp)))


Sigma <- Reduce(`+`, lapply(aux_list, solve))


True_Sigma <- gets_true_cov_mat(graph = graph_initial, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)


p <- graph_initial$plot_function(True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
graph_initial$plot_function(Sigma[,2], p = p, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)


