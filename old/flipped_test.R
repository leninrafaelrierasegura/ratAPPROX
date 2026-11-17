# this is one is ok the next is to see what Alex says

library(Matrix)
library(MetricGraph)
library(viridis)
library(dplyr)

# get the functions from the Rmd files
knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R"))
knitr::purl(here::here("matern_functions.Rmd"), output = here::here("matern_functions.R"))

# load the functions
source(here::here("all_functions.R"))
source(here::here("matern_functions.R"))

# parameters
h <- 0.05
kappa <- 10
sigma <- 0.8
alpha <- 2.8
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000

type_rational_approx = "brasil"
type_interp = "spline"

FLIPPED <- FALSE


# build a graph with a mesh
graph_initial <- gets.graph.tadpole(flip_edge = FLIPPED)
graph <- graph_initial$clone()


graph_initial$build_mesh(h = h)
graph_initial$compute_fem()

# mesh_locations <- graph_initial$get_mesh_locations() %>% 
#   as.data.frame() %>% 
#   mutate(dummy = 1) %>%
#   rename(edge_number = V1, distance_on_edge = V2)

mesh_locations <- graph_initial$get_mesh_locations() %>%
  as.data.frame() %>%
  mutate(dummy = 1, mesh_row = row_number()) %>%
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(
  data = mesh_locations,
  normalized = TRUE)

# add observations as vertices
graph$observation_to_vertex()


# RAT_Sigma <- rat_covariance(
#   graph = graph, 
#   kappa = kappa, 
#   tau = tau, 
#   alpha = alpha, 
#   m = m,
#   type_rational_approx = type_rational_approx,
#   type_interp = type_interp,
#   build_cov = TRUE)

obs <- as.data.frame(graph$get_data())
mesh_order <- order(obs$mesh_row)

# follows the observation order returned by graph$get_data()
AQ <- rat_covariance(
  graph = graph,
  kappa = kappa,
  tau = tau,
  alpha = alpha,
  m = m,
  type_rational_approx = type_rational_approx,
  type_interp = type_interp,
  build_cov = FALSE
)

A_mesh <- AQ$A[mesh_order, , drop = FALSE]

RAT_Sigma <- A_mesh %*% solve(AQ$Q, t(A_mesh))


op <- matern.operators(alpha = alpha, 
                      kappa = kappa, 
                      tau = tau,
                      m = m, 
                      graph = graph_initial)

# FEM_Sigma follows graph_initial$mesh$V
FEM_Sigma <- covariance_mesh(op)

True_Sigma <- gets_true_cov_mat(
  graph = graph_initial, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  n.overkill = n.overkill)

q <- graph_initial$plot_function(X = FEM_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "FEM", showlegend = TRUE)
graph_initial$plot_function(X = RAT_Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "RAT", showlegend = TRUE)



q <- graph_initial$plot_function(X = diag(FEM_Sigma), type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "FEM", showlegend = TRUE)
graph_initial$plot_function(X = diag(RAT_Sigma), p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "RAT", showlegend = TRUE)



L_2_error <- sqrt(as.double(t(graph_initial$mesh$weights)%*%(FEM_Sigma - RAT_Sigma)^2%*%graph_initial$mesh$weights))
L_2_error
L_2_error <- sqrt(as.double(t(graph_initial$mesh$weights)%*%(True_Sigma - RAT_Sigma)^2%*%graph_initial$mesh$weights))
L_2_error
L_2_error <- sqrt(as.double(t(graph_initial$mesh$weights)%*%(FEM_Sigma - True_Sigma)^2%*%graph_initial$mesh$weights))
L_2_error
