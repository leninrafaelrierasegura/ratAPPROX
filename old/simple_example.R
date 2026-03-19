
library(MetricGraph)
library(Matrix)
library(dplyr)


source("old/functions_for_simple_example.R")

# parameters
h <- 0.2
kappa <- 5
sigma <- 0.8
alpha <- 1.001
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000


FLIPPED <- TRUE


# build a graph with a mesh
graph_initial <- gets.graph.tadpole(flip_edge = FLIPPED)
graph <- graph_initial$clone()




graph_initial$build_mesh(h = h)
graph_initial$compute_fem()


# add mesh locations as observations
graph$add_observations(
  data = graph_initial$get_mesh_locations() %>% 
    as.data.frame() %>% 
    mutate(y = 1) %>% 
    rename(edge_number = V1, distance_on_edge = V2) %>% 
    mutate(dummy = 1:length(y)),
  edge_number = "edge_number",
  distance_on_edge = "distance_on_edge",
  data_coords = "PtE",
  normalized = TRUE, 
  clear_obs = TRUE)

# add observations as vertices
graph$observation_to_vertex()


Approx_Sigma <- gets_cov_mat_rat_approx_alpha_1_to_2(
  graph = graph, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  m = m)


graph_true <- gets.graph.tadpole(flip_edge = FALSE)
graph_true$build_mesh(h = h)
graph_true$compute_fem()
True_Sigma <- gets_true_cov_mat(graph = graph_true, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)

# mesh points in tail edge
nt <- sum(graph_true$mesh$VtE[,1] == 1) + 1

my_order <- if (FLIPPED) c(nt, 1, (nt-1):2, (nt+1):graph$nV) else 1:graph$nV
Approx_Sigma <- Approx_Sigma[my_order, my_order]



q <- graph_true$plot_function(X = True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)

graph_true$plot_function(X = Approx_Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)



L_2_error = sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graph_true$mesh$weights))
print(L_2_error)













