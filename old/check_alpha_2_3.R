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
h <- 0.001
kappa <- 24
sigma <- 1
alpha <- 2.5
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
n.overkill <- 1000


# build a graph with a mesh
graphFEM <- gets.graph.tadpole(flip_edge = FALSE)
graph <- graphFEM$clone()




graphFEM$build_mesh(h = h)

mesh_XY <- graphFEM$mesh$V %>% as.data.frame()  

mesh_locations <- graphFEM$get_mesh_locations() %>% 
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


Approx_Sigma <- getsCovarianceMatrixForRationalApproximationForAlphaBetweenTwoAndThree(
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


#Approx_Sigma <- Approx_Sigma[idx, idx]


True_Sigma <- gets_true_cov_mat(graph = graphFEM, 
                                kappa = kappa, 
                                tau = tau, 
                                alpha = alpha, 
                                n.overkill = n.overkill)


op = matern.operators(alpha = alpha, 
                      kappa = kappa, 
                      tau = tau,
                      m = m, 
                      graph = graphFEM)

appr_cov_mat = covariance_mesh(op)


# q <- graphFEM$plot_function(X = True_Sigma[,1], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)
# 
# graphFEM$plot_function(X = Approx_Sigma[,1], p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
#   graphFEM$plot_function(X = appr_cov_mat[,1], p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)
# 

q <- graphFEM$plot_function(X = diag(True_Sigma), type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)

graphFEM$plot_function(X = diag(Approx_Sigma), p = q, type = "plotly", line_color = "blue", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE) %>%
  graphFEM$plot_function(X = diag(appr_cov_mat), p = ., type = "plotly", line_color = "green", interpolate_plot = FALSE, name = "Approx", showlegend = TRUE)



L_2_error = sqrt(as.double(t(graphFEM$mesh$weights)%*%(True_Sigma - Approx_Sigma)^2%*%graphFEM$mesh$weights))
print(L_2_error)

