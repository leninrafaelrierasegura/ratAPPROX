# this is one is ok the next is to see what Alex says

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
h <- 0.1
kappa <- 1
sigma <- 1
alpha <- 0.9
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000

# This is a function from MetricGraph package that returns a list of edges
edges <- logo_lines() 
# Create a new graph object
graph_initial <- metric_graph$new(edges = edges, perform_merges = TRUE) 
# Build the mesh

graph <- graph_initial$clone()


graph_initial$build_mesh(h = h)


mesh_locations <- graph_initial$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(dummy = 1) %>%
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(
  data = mesh_locations,
  normalized = TRUE)

# add observations as vertices
graph$observation_to_vertex()


data <- graph$get_data()

Approx_Sigma <- rat_covariance(
  graph = graph, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  m = m)



# add mesh locations as observations
graph$add_observations(
  data = data %>% mutate(cov = Approx_Sigma[,100]),
  normalized = TRUE, clear_obs = TRUE)

# # add observations as vertices
# graph$observation_to_vertex()


graph$plot_function(data = "cov", 
                    type = "plotly", 
                    interpolate_plot = FALSE, 
                    vertex_size = 0)


op <- matern.operators(alpha = alpha, 
                      kappa = kappa, 
                      tau = tau,
                      m = m, 
                      graph = graph_initial)

FEM_Sigma <- covariance_mesh(op)


graph_initial$plot_function(X = Approx_Sigma[,100],  
                            type = "plotly", 
                            line_color = "green", 
                            interpolate_plot = FALSE, 
                            name = "FEM", 
                            showlegend = TRUE)





