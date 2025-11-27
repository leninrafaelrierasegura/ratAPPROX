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

graph_initial$build_mesh(h = 0.6)
# to get the weights
graph_initial$compute_fem()
# add mesh locations as observations
graph$add_observations(data = graph_initial$get_mesh_locations() %>% 
                         as.data.frame() %>% 
                         mutate(y = 1) %>% 
                         rename(edge_number = V1, distance_on_edge = V2),
                       edge_number = "edge_number",
                       distance_on_edge = "distance_on_edge",
                       data_coords = "PtE",
                       normalized = TRUE, 
                       clear_obs = TRUE)
# add observations as vertices
graph$observation_to_vertex()

# parameters
kappa <- 10
sigma <- 1
alpha <- 1.9
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 100

Sigma <- gets_cov_mat_rat_approx_alpha_1_to_2(graph = graph, 
                                              kappa = kappa, 
                                              tau = tau, 
                                              alpha = alpha, 
                                              m = m)

True_Sigma <- gets_true_cov_mat(graph = graph_initial, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)

q <- graph_initial$plot_function(True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
graph_initial$plot_function(Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)

L_2_error = sqrt(as.double(t(graph_initial$mesh$weights)%*%(True_Sigma - Sigma)^2%*%graph_initial$mesh$weights))
print(L_2_error)


  