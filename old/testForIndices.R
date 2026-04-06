capture.output(
  knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R")),
  file = here::here("old/purl_log.txt")
)
capture.output(
  knitr::purl(here::here("matern_functions.Rmd"), output = here::here("matern_functions.R")),
  file = here::here("old/purl_log.txt")
)

source(here::here("all_functions.R"))
source(here::here("matern_functions.R"))

theta1 <- seq(from=pi,to=0,length.out = 1000)
theta2 <- seq(from=-pi,to=0,length.out = 1000)

edge1 <- cbind(1/pi+cos(theta1)/pi,sin(theta1)/pi)
edge2 <- rbind(c(0,0),c(2/pi,0))
edge3 <- cbind(1/pi+cos(theta2)/pi,sin(theta2)/pi)
edges <- list(edge1, edge2, edge3)
graph <- metric_graph$new(edges = edges)
graph$set_manual_edge_lengths(edge_lengths = c(1,2/pi,1))

kappa <- 4
sigma <- 1
alpha <- 1
m <- 4
nu <- alpha - 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))


Q_tilde <- computesListOfMatricesQTildeUnconstraint(p = 0,
                                                     kappa = kappa, 
                                                     alpha = alpha, 
                                                     edgeLengths = graph$edge_lengths)[[1]]

graph$buildC()
C <- graph$C[c(2,3,5,6), c(1,3,5,7,9,11)]

K <- buildKirchooffConditioningMatrixCaseAlphaEqualOne(graph)
CoB <- MetricGraph:::c_basis2(K)

CoB$T <- t(CoB$T)

Tc <- CoB$T[-c(1:length(CoB$S)), ]

graph$buildC()
graph$C
graph$CoB


Q_tilde_star <- Tc %*% Q_tilde %*% t(Tc)
vertex_loc <- graph$V %>% 
  as.data.frame() %>%
  rename(coord_x = X, coord_y = Y) %>%
  mutate(dummy = 1)

# add mesh locations as observations
graph$add_observations(
  data = vertex_loc,
  coord_x = "coord_x",
  coord_y = "coord_y",
  data_coords = "spatial",
)

# add observations as vertices
graph$observation_to_vertex()

index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
index.obs_i <- (index.obs_i - 1)/2 + 1

A <- t(Tc)[index.obs_i, ]

Q_tilde_star_UU <- Q_tilde_star

Sigma_build <- A %*%  solve(Q_tilde_star_UU, t(A))





Approx_Sigma <- rat_covariance(
  graph = graph, 
  kappa = kappa, 
  tau = tau, 
  alpha = alpha, 
  m = m)


Approx_Sigma

Sigma_build


sum(Approx_Sigma-Sigma_build)