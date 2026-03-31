
set.seed(1982)
#devtools::load_all()
library(rSPDE)

edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edges = list(edge1, edge2, edge3, edge4)
graph <- metric_graph$new(edges = edges)

graph_clone <- graph$clone()

graph$build_mesh(h = 0.01)

alpha <- 1
nu <- alpha - 1/2
sigma <- 1.3
range <- 0.15
kappa <- sqrt(8*nu)/range
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))


rspde.order <- 2



op <- matern.operators(alpha = alpha,
                       kappa = kappa,
                       tau = tau,
                       parameterization = "spde",
                       m = rspde.order,
                       graph = graph)


u <- simulate(op)

obs.per.edge <- 200
obs.loc <- NULL
for(i in 1:graph$nE) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge)))
}
n.obs <- obs.per.edge*graph$nE
A <- graph$fem_basis(obs.loc)

sigma_e <- 0.5
Y <- as.vector(A %*% u + sigma_e * rnorm(n.obs))
df_data <- data.frame(y = Y, edge_number = obs.loc[,1],
                      distance_on_edge = obs.loc[,2])
graph_clone$add_observations(data = df_data, normalized = TRUE)

#undebug(graph_lme)



fit <- graph_lme(y ~ -1, graph=graph_clone, BC = 0,
                 model = list(type = "WhittleMatern", alpha = 1))

summary(fit)



source("lenin/minimal_functions.R")

#undebug(graph_lme_minimal_v2)


fit <- graph_lme_minimal_v2(y ~ -1, graph=graph_clone, BC = 0,
                         model = list(type = "WhittleMatern", alpha = 1, version = 2))

summary(fit)



fit <- graph_lme_minimal(y ~ -1, graph=graph_clone, BC = 0,
                 model = list(type = "WhittleMatern", alpha = 1, version = 1))

summary(fit)

precomp_data <- fit$precomp_data

# using true parameters
eval_likelihood_alpha1(
  sigma_e     = sigma_e,
  tau         = tau,
  kappa       = kappa,
  Y           = Y,
  graph       = graph_clone,
  precomp_data = precomp_data,
  BC          = 0
)

# using fitted parameters
eval_likelihood_alpha1(
  sigma_e      = fit$coeff$measurement_error,
  tau          = fit$coeff$random_effects[["tau"]],
  kappa        = fit$coeff$random_effects[["kappa"]],
  graph        = graph_clone,
  precomp_data = precomp_data,
  BC           = 0
)



graph_clone$observation_to_vertex()

data <- graph_clone$get_data()

Y <- data$y

sigma_e      = fit$coeff$measurement_error
tau          = fit$coeff$random_effects[["tau"]]
kappa        = fit$coeff$random_effects[["kappa"]]

theta <- c(log(sigma_e), log(1/tau), log(kappa))

eval_likelihood_alpha1_v2(theta = theta,
                          graph = graph_clone,
                          X_cov = matrix(0, ncol = 0, nrow = 1),
                          y = Y,
                          repl = NULL,
                          BC = 0,
                          parameterization = "spde")







