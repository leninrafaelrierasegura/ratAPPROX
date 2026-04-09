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
h <- 0.01
kappa <- 5
sigma <- 0.8
alpha <- 3
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000


FLIPPED <- FALSE


# build a graph with a mesh
graph_initial <- gets.graph.tadpole(flip_edge = FLIPPED)
graph <- graph_initial$clone()


graph_initial$build_mesh(h = h)

mesh_XY <- graph_initial$mesh$V %>% as.data.frame()

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

order <- data %>%
  as.data.frame() %>%
  select(.coord_x, .coord_y) %>%
  rename(X = .coord_x, Y = .coord_y)



Approx_Sigma <- rat_covariance(
  graph = graph, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  m = m,
  type_rational_approx = "chebfun",
  type_interp = "spline",
  build_cov = TRUE)

digits <- 10

idx <- match(
  paste(round(mesh_XY$X, digits), round(mesh_XY$Y, digits)),
  paste(round(order$X, digits), round(order$Y, digits))
)


Approx_Sigma <- Approx_Sigma[idx, idx]

graph_true <- gets.graph.tadpole(flip_edge = FALSE)
graph_true$build_mesh(h = h)
graph_true$compute_fem()

True_Sigma <- gets_true_cov_mat(
  graph = graph_true, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  n.overkill = n.overkill)

# mesh points in tail edge
# nt <- sum(graph_true$mesh$VtE[,1] == 1) + 1
# 
# if (!FLIPPED) {my_order <- 1:graph$nV
# } else if (alpha > 0.5 && alpha <= 1) {my_order <- c(2, 1, nt:3, (nt+1):graph$nV)
#   } else if (alpha > 1 && alpha <= 2) {my_order <- c(nt, 1, (nt-1):2, (nt+1):graph$nV)
#     } else {stop("alpha outside supported range")
#       }


#Approx_Sigma <- Approx_Sigma[my_order, my_order]


op <- matern.operators(alpha = alpha, 
                      kappa = kappa, 
                      tau = tau,
                      m = m, 
                      graph = graph_true)

FEM_Sigma <- covariance_mesh(op)

q <- graph_true$plot_function(X = True_Sigma[,2], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
graph_true$plot_function(X = Approx_Sigma[,2], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
  graph_true$plot_function(X = FEM_Sigma[,2], p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "FEM", showlegend = TRUE)



# q <- graph_true$plot_function(X = diag(True_Sigma), type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
# graph_true$plot_function(X = diag(Approx_Sigma), p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
#   graph_true$plot_function(X = diag(FEM_Sigma), p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "FEM", showlegend = TRUE)
# 


L_2_error_RAT <- sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graph_true$mesh$weights))
L_2_error_FEM <- sqrt(as.double(t(graph_true$mesh$weights)%*%(True_Sigma - FEM_Sigma)^2%*%graph_true$mesh$weights))


cat(sprintf("L2 error (RAT-TRUE): %.12f\nL2 error (FEM-TRUE): %.12f\n", L_2_error_RAT, L_2_error_FEM))
