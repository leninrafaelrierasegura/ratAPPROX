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
h <- 0.2
kappa <- 5
sigma <- 0.8
alpha <- 1
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000


FLIPPED <- FALSE


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

In <- Diagonal(graph$nV)


Approx_Sigma <- solve(MetricGraph:::Qalpha1(theta = c(tau, kappa),
                        graph = graph,
                        BC = 300,
                        build = TRUE), In)




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

my_order <- if (FLIPPED) c(2, 1, nt:3, (nt+1):graph$nV) else 1:graph$nV
Approx_Sigma <- Approx_Sigma[my_order, my_order]


# order_tmp <- if (FLIPPED) c(2, 6, 5, 4, 3, 1, 7:graph$nV) else 1:graph$nV
# True_Sigma <- True_Sigma[order_tmp,order_tmp]


op = matern.operators(alpha = alpha, 
                      kappa = kappa, 
                      tau = tau,
                      m = m, 
                      graph = graph_true)

appr_cov_mat = covariance_mesh(op)


q <- graph_true$plot_function(X = True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)

graph_true$plot_function(X = Approx_Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
  graph_true$plot_function(X = appr_cov_mat[,2], p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)


q <- graph_true$plot_function(X = diag(True_Sigma), type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)

graph_true$plot_function(X = diag(Approx_Sigma), p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
  graph_true$plot_function(X = diag(appr_cov_mat), p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)



L_2_error = sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graph_true$mesh$weights))
print(L_2_error)
