library(Matrix)
library(MetricGraph)
library(viridis)

# get the functions from the Rmd files
knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R"))
knitr::purl(here::here("matern_functions.Rmd"), output = here::here("matern_functions.R"))

# load the functions
source(here::here("all_functions.R"))
source(here::here("matern_functions.R"))

# parameters
h <- 0.2
kappa <- 5
sigma <- 0.8
alpha <- 1.4
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000

# build a graph with a mesh
graph_initial <- gets.graph.tadpole(flip_edge = TRUE)
graph <- graph_initial$clone()

graph_initial$build_mesh(h = h)
# to get the weights
graph_initial$compute_fem()
# add mesh locations as observations
graph$add_observations(
  data = graph_initial$get_mesh_locations() %>% 
    as.data.frame() %>% 
    mutate(y = 1) %>% 
    rename(edge_number = V1, distance_on_edge = V2)%>%mutate(dummy=1:length(y)),
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



order_tmp <- graph$get_data()$dummy
Approx_Sigma[order_tmp,] <- Approx_Sigma[1:nrow(Approx_Sigma),]

graph_true <- gets.graph.tadpole(flip_edge = FALSE)
graph_true$build_mesh(h = h)
graph_true$compute_fem()
True_Sigma <- gets_true_cov_mat(graph = graph_true, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)

q <- graph_true$plot_function(X = True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)

graph_initial$plot_function(X = Approx_Sigma[,1], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)


Approx_Sigma[order_tmp, ] <- Approx_Sigma[,order_tmp]


L_2_error = sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graph_true$mesh$weights))
print(L_2_error)


# graph_true <- gets.graph.tadpole(flip_edge = FALSE)
# graph_true$build_mesh(h = h)
# VtE<- graph_true$mesh$VtE
# ord <- c(2,1,max(which(VtE[,1]==1)):3, (max(which(VtE[,1]==1))+1):nrow(VtE))
# Approx_Sigma <- Approx_Sigma[ord, ord]
# # cbind(graph$V[ord,], graph_true$mesh$V)
# # cbind(graph_initial$mesh$VtE, graph_true$mesh$VtE)
# graph_true$compute_fem()
# True_Sigma <- gets_true_cov_mat(graph = graph_true, 
#                                 kappa = kappa, 
#                                 tau = tau, 
#                                 alpha = alpha, 
#                                 n.overkill = n.overkill)
# 
# q <- graph_true$plot_function(X = True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
# graph_true$plot_function(X = Approx_Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)
# 
# L_2_error = sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graph_true$mesh$weights))
# print(L_2_error)

